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
    gzip_content_fingerprint,
    large_file_fingerprint,
    safe_remove_owned,
    small_file_fingerprint,
	tool_identity,
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


def test_path_fingerprints_are_stable_across_equivalent_directory_aliases(tmp_path):
    physical = tmp_path / "physical"
    physical.mkdir()
    alias = tmp_path / "alias"
    alias.symlink_to(physical, target_is_directory=True)
    (physical / "input.txt").write_text("content\n")

    physical_fingerprint = small_file_fingerprint(physical / "input.txt")
    alias_fingerprint = small_file_fingerprint(alias / "input.txt")

    assert alias_fingerprint == physical_fingerprint
    assert CacheManager(str(alias / "run")).cache_root == str(physical / "run.cache")
    assert OutputRequirement("out", str(alias / "out.txt")).normalized_path() == str(
        physical / "out.txt"
    )


def test_gzip_content_fingerprint_ignores_mtime_but_detects_content(tmp_path):
    import gzip

    path = tmp_path / "reads.fq.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("@read\nACGT\n+\nIIII\n")
    first = gzip_content_fingerprint(path)

    stat = path.stat()
    os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000_000))
    assert gzip_content_fingerprint(path) == first

    with gzip.open(path, "wt") as handle:
        handle.write("@read\nTGCA\n+\nIIII\n")
    assert gzip_content_fingerprint(path)["crc32"] != first["crc32"]


def test_tool_identity_is_resolved_once_per_process(monkeypatch):
    calls = []
    tool_identity.cache_clear()
    monkeypatch.setattr("CRISPRSCope.cache.shutil.which", lambda _name: "/tools/example")
    monkeypatch.setattr(
        "CRISPRSCope.cache.subprocess.check_output",
        lambda *args, **kwargs: calls.append(args[0]) or "example 1.0  \n",
    )
    try:
        first = tool_identity(("example-cache-test", "--version"))
        second = tool_identity(("example-cache-test", "--version"))
    finally:
        tool_identity.cache_clear()
    assert first == second == {"path": "/tools/example", "version": "example 1.0"}
    assert calls == [("example-cache-test", "--version")]


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
    disabled._read_record = lambda _path: (_ for _ in ()).throw(
        AssertionError("disabled mode must not read cache records")
    )
    assert disabled.load("stage") is None
    disabled.commit(disabled_record, [requirement])
    assert not os.path.exists(disabled.record_path("disabled"))


def test_malformed_cell_counts_and_gzip_outputs_invalidate_cleanly(tmp_path):
    manager = CacheManager(str(tmp_path / "run"))
    counts = tmp_path / "counts.txt"
    counts.write_text("cellA\t1\n")
    fastq = tmp_path / "reads.fq.gz"
    import gzip
    with gzip.open(fastq, "wt") as handle:
        handle.write("")
    requirements = (
        OutputRequirement(
            "counts", str(counts), strategy="sha256", validator="cell_counts"
        ),
        OutputRequirement(
            "fastq", str(fastq), strategy="stat", allow_empty=True,
            validator="gzip",
        ),
    )
    record = manager.new_record("validated", algorithm_version=1)
    manager.commit(record, requirements)

    counts.write_text("cellA\tnot-an-integer\n")
    decision = manager.evaluate(record, requirements)
    assert decision.status == "invalid"
    assert decision.reasons == ("output_invalid:counts",)

    counts.write_text("cellA\t1\n")
    fastq.write_text("not gzip\n")
    decision = manager.evaluate(record, requirements)
    assert decision.status == "invalid"
    assert decision.reasons == ("output_invalid_gzip:fastq",)


def test_fastq_validator_accepts_gzip_or_plain_crispresso_output(tmp_path):
    manager = CacheManager(str(tmp_path / "run"))
    plain = tmp_path / "plain.fastq.gz"
    plain.write_text("@plain\nACGT\n+\nIIII\n")
    compressed = tmp_path / "compressed.fastq.gz"
    import gzip
    with gzip.open(compressed, "wt") as handle:
        handle.write("@gzip\nACGT\n+\nIIII\n")
    requirements = (
        OutputRequirement("plain", str(plain), strategy="stat", validator="fastq"),
        OutputRequirement(
            "compressed", str(compressed), strategy="stat", validator="fastq"
        ),
    )
    record = manager.new_record("crispresso_fastq", algorithm_version=1)
    manager.commit(record, requirements)
    assert manager.evaluate(record, requirements).is_hit

    plain.write_text("not a FASTQ\n")
    decision = manager.evaluate(record, requirements)
    assert decision.status == "invalid"
    assert decision.reasons == ("output_invalid_fastq:plain",)


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


def test_safe_remove_owned_accepts_equivalent_directory_aliases(tmp_path):
    physical_root = tmp_path / "physical-root"
    physical_root.mkdir()
    alias_root = tmp_path / "alias-root"
    alias_root.symlink_to(physical_root, target_is_directory=True)
    candidate = alias_root / "expected.txt"
    candidate.write_text("remove\n")

    assert safe_remove_owned(
        candidate,
        allowed_root=physical_root,
        expected_name="expected.txt",
    )
    assert not (physical_root / "expected.txt").exists()


def test_trusted_stage_inventory_ignores_corrupt_and_misnamed_records(tmp_path, caplog):
    manager = CacheManager(str(tmp_path / "run"))
    trusted = manager.new_record("stage", "ampA", algorithm_version=1)
    manager.commit(trusted, ())
    corrupt = tmp_path / "run.cache" / "stage" / "corrupt.json"
    corrupt.write_text("{")

    records = manager.trusted_stage_records("stage")

    assert [record.scope for record in records] == ["ampA"]
    assert "Ignoring untrusted cache record" in caplog.text
    assert manager.remove_record(records[0])
    assert not os.path.exists(manager.record_path("stage", "ampA"))
