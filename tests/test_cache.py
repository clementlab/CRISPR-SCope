import json
import os

import pytest

from CRISPRSCope.cache import (
    CacheConfig,
    CacheManager,
    CacheMode,
    OutputRequirement,
    OutputRootLock,
    canonical_digest,
    large_file_fingerprint,
    safe_remove_owned,
    small_file_fingerprint,
)


def test_canonical_digest_normalizes_mapping_and_set_order():
    left = {"b": {"cellB", "cellA"}, "a": [2, 1]}
    right = {"a": [2, 1], "b": {"cellA", "cellB"}}
    assert canonical_digest(left) == canonical_digest(right)
    assert canonical_digest(left) != canonical_digest({"a": [1, 2], "b": {"cellA", "cellB"}})


def test_cache_config_validates_modes():
    assert CacheConfig.from_value("AUTO").mode is CacheMode.AUTO
    assert CacheConfig.from_value(" refresh ").mode is CacheMode.REFRESH
    assert CacheConfig.from_value("disabled").mode is CacheMode.DISABLED
    with pytest.raises(ValueError, match="Invalid value for cache_mode"):
        CacheConfig.from_value("sometimes")


def test_small_and_large_fingerprints_use_declared_strategies(tmp_path, monkeypatch):
    path = tmp_path / "input.fastq"
    path.write_text("ACGT\n")
    assert small_file_fingerprint(path)["sha256"]

    def forbidden_open(*_args, **_kwargs):
        raise AssertionError("large-file fingerprint must not read the payload")

    monkeypatch.setattr("builtins.open", forbidden_open)
    fingerprint = large_file_fingerprint(path)
    assert fingerprint["strategy"] == "stat"
    assert fingerprint["size"] == 5
    assert "sha256" not in fingerprint


def test_cache_record_round_trip_and_output_invalidation(tmp_path):
    output_root = str(tmp_path / "run")
    output = tmp_path / "result.tsv"
    output.write_text("name\tvalue\na\t1\n")
    manager = CacheManager(output_root)
    record = manager.new_record(
        "example",
        algorithm_version=1,
        inputs={"control": "A"},
        parameters={"threshold": 2},
    )
    requirement = OutputRequirement(
        "table",
        str(output),
        strategy="sha256",
        validator="tsv",
        required_header=("name", "value"),
    )

    assert manager.evaluate(record, [requirement]).status == "miss"
    committed = manager.commit(record, [requirement])
    decision = manager.evaluate(
        manager.new_record(
            "example",
            algorithm_version=1,
            inputs={"control": "A"},
            parameters={"threshold": 2},
        ),
        [requirement],
    )
    assert decision.is_hit
    assert decision.record.cache_key == committed.cache_key

    output.write_text("wrong\theader\na\t1\n")
    invalid = manager.evaluate(record, [requirement])
    assert invalid.status == "invalid"
    assert invalid.reasons == ("output_invalid_header:table",)


def test_changed_inputs_and_unsupported_or_truncated_records_are_invalid(tmp_path):
    manager = CacheManager(str(tmp_path / "run"))
    output = tmp_path / "output.txt"
    output.write_text("ok\n")
    requirement = OutputRequirement("output", str(output), strategy="sha256")
    original = manager.new_record("stage", algorithm_version=1, inputs={"value": 1})
    manager.commit(original, [requirement])

    changed = manager.new_record("stage", algorithm_version=1, inputs={"value": 2})
    assert manager.evaluate(changed, [requirement]).reasons == ("cache_key_changed",)

    record_path = manager.record_path("stage")
    payload = json.loads(open(record_path, encoding="utf-8").read())
    payload["schema_version"] = 999
    with open(record_path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle)
    assert manager.evaluate(original, [requirement]).status == "invalid"

    with open(record_path, "w", encoding="utf-8") as handle:
        handle.write("{")
    assert manager.evaluate(original, [requirement]).status == "invalid"


def test_refresh_and_disabled_modes_do_not_hit(tmp_path):
    output = tmp_path / "output.txt"
    output.write_text("ok\n")
    requirement = OutputRequirement("output", str(output), strategy="sha256")
    auto = CacheManager(str(tmp_path / "run"))
    record = auto.new_record("stage", algorithm_version=1)
    auto.commit(record, [requirement])

    refresh = CacheManager(str(tmp_path / "run"), CacheConfig(CacheMode.REFRESH))
    assert refresh.evaluate(refresh.new_record("stage", algorithm_version=1), [requirement]).status == "refresh"

    disabled = CacheManager(str(tmp_path / "run"), CacheConfig(CacheMode.DISABLED))
    disabled_record = disabled.new_record("disabled", algorithm_version=1)
    assert disabled.evaluate(disabled_record, [requirement]).status == "disabled"
    disabled.commit(disabled_record, [requirement])
    assert not os.path.exists(disabled.record_path("disabled"))


def test_per_scope_record_paths_do_not_collide(tmp_path):
    manager = CacheManager(str(tmp_path / "run"))
    assert manager.record_path("stage", "amp/A") != manager.record_path("stage", "amp_A")
    assert manager.record_path("stage", "amp/A").endswith(".json")


def test_output_root_lock_fails_immediately_for_second_owner(tmp_path):
    first = OutputRootLock(str(tmp_path / "run")).acquire()
    try:
        with pytest.raises(RuntimeError, match="Another CRISPRSCope process"):
            OutputRootLock(str(tmp_path / "run")).acquire()
    finally:
        first.release()

    with OutputRootLock(str(tmp_path / "run")):
        pass


def test_safe_remove_owned_rejects_outside_paths_and_external_symlinks(tmp_path):
    owned = tmp_path / "owned"
    owned.mkdir()
    inside = owned / "expected.txt"
    inside.write_text("remove\n")
    assert safe_remove_owned(inside, allowed_root=owned, expected_name="expected.txt")
    assert not inside.exists()

    outside = tmp_path / "outside.txt"
    outside.write_text("keep\n")
    with pytest.raises(ValueError, match="outside cache-owned root"):
        safe_remove_owned(outside, allowed_root=owned)
    assert outside.exists()

    link = owned / "external-link"
    link.symlink_to(outside)
    with pytest.raises(ValueError, match="symlink outside"):
        safe_remove_owned(link, allowed_root=owned)
    assert outside.exists()
