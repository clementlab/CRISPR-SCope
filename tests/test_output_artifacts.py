import json

import pytest

from CRISPRSCope.output_artifacts import OutputContext, OutputManifest, OutputSpec


def test_context_preserves_existing_paths_and_expands_plot_formats(tmp_path):
    context = OutputContext(str(tmp_path / "run"), h5ad_output=str(tmp_path / "custom.h5ad"))

    assert context.path("editing_summary") == str(tmp_path / "run.editingSummary.txt")
    assert context.path("h5ad") == str(tmp_path / "custom.h5ad")
    assert context.paths("editing_rate_unconditional_permutation_plot") == (
        str(tmp_path / "run.12_EditingRateUnconditionalPermutation.png"),
        str(tmp_path / "run.12_EditingRateUnconditionalPermutation.pdf"),
    )
    assert context.data_links("editing_rate_unconditional_permutation_plot") == [
        ("Unconditional permutation summary", str(tmp_path / "run.editingRateUnconditionalPermutation.txt")),
        ("Unconditional permutation simulations", str(tmp_path / "run.editingRateUnconditionalPermutationSimulations.txt")),
    ]


def test_context_rejects_unknown_and_duplicate_keys(tmp_path):
    context = OutputContext(str(tmp_path / "run"))
    with pytest.raises(KeyError, match="Unknown output artifact key"):
        context.path("missing")

    duplicate_specs = (
        OutputSpec("duplicate", ".one", "table"),
        OutputSpec("duplicate", ".two", "table"),
    )
    with pytest.raises(RuntimeError, match="must be unique"):
        OutputContext(str(tmp_path / "run"), specs=duplicate_specs)


def test_intermediate_family_delegates_to_existing_filename_builder(tmp_path):
    context = OutputContext(str(tmp_path / "run"))
    calls = []

    def filename_builder(**kwargs):
        calls.append(kwargs)
        return "expected-path"

    assert context.family("amplicon_fastq").resolve(
        filename_builder, stage=3, tag="reads_all_cells", amplicon="ampA", read="r1"
    ) == "expected-path"
    assert calls == [{"stage": 3, "tag": "reads_all_cells", "amplicon": "ampA", "read": "r1"}]


def test_remove_only_declared_artifacts(tmp_path):
    context = OutputContext(str(tmp_path / "run"))
    declared = tmp_path / "run.14_EditingRateDepthStability.png"
    declared.write_text("stale")
    unrelated = tmp_path / "run.keep-me.txt"
    unrelated.write_text("keep")

    assert context.remove(("editing_rate_depth_stability_plot",)) == [str(declared)]
    assert not declared.exists()
    assert unrelated.exists()


def test_manifest_records_completed_and_failed_runs_atomically(tmp_path):
    context = OutputContext(str(tmp_path / "run"))
    (tmp_path / "run.editingSummary.txt").write_text("summary")
    manifest = OutputManifest(context)
    manifest.mark_written("editing_summary")
    manifest.mark_skipped("editing_rate_depth_stability_plot", "analysis is disabled")
    manifest.complete()
    manifest_path = manifest.write()

    payload = json.loads(open(manifest_path, encoding="utf-8").read())
    assert payload["schema_version"] == 1
    assert payload["status"] == "completed"
    assert [artifact["key"] for artifact in payload["artifacts"]][:2] == [
        "editing_summary",
        "filtered_editing_summary",
    ]
    assert next(item for item in payload["artifacts"] if item["key"] == "output_manifest")["status"] == "written"
    skipped = next(item for item in payload["artifacts"] if item["key"] == "editing_rate_depth_stability_plot")
    assert skipped["status"] == "skipped"
    assert skipped["reason"] == "analysis is disabled"
    assert not list(tmp_path.glob(".run.outputManifest.json.*.tmp"))

    failed = OutputManifest(context)
    failed.set_stage("run_crispresso")
    error = ValueError("simulated failure")
    failed.fail("run_crispresso", error)
    failure_payload = failed.as_dict()
    assert failure_payload["status"] == "failed"
    assert failure_payload["failure"] == {
        "stage": "run_crispresso",
        "type": "ValueError",
        "message": "simulated failure",
    }


def test_cli_main_writes_requested_manifest_on_success_and_preserves_failure(tmp_path, monkeypatch):
    from CRISPRSCope import cli

    context = OutputContext(str(tmp_path / "success"))

    def successful_pipeline():
        cli._ACTIVE_OUTPUT_MANIFEST = OutputManifest(context)
        cli._ACTIVE_OUTPUT_MANIFEST.set_stage("write_summary")
        (tmp_path / "success.editingSummary.txt").write_text("summary")

    monkeypatch.setattr(cli, "_main_impl", successful_pipeline)
    assert cli.main() is None
    success_payload = json.loads((tmp_path / "success.outputManifest.json").read_text())
    assert success_payload["status"] == "completed"
    assert next(item for item in success_payload["artifacts"] if item["key"] == "editing_summary")["status"] == "written"

    failure_context = OutputContext(str(tmp_path / "failure"))

    def failing_pipeline():
        cli._ACTIVE_OUTPUT_MANIFEST = OutputManifest(failure_context)
        cli._ACTIVE_OUTPUT_MANIFEST.set_stage("run_crispresso")
        raise RuntimeError("expected pipeline failure")

    monkeypatch.setattr(cli, "_main_impl", failing_pipeline)
    with pytest.raises(RuntimeError, match="expected pipeline failure"):
        cli.main()
    failure_payload = json.loads((tmp_path / "failure.outputManifest.json").read_text())
    assert failure_payload["status"] == "failed"
    assert failure_payload["failure"] == {
        "stage": "run_crispresso",
        "type": "RuntimeError",
        "message": "expected pipeline failure",
    }


def test_cli_main_does_not_write_manifest_without_an_active_manifest(tmp_path, monkeypatch):
    from CRISPRSCope import cli

    monkeypatch.setattr(cli, "_main_impl", lambda: None)
    assert cli.main() is None
    assert not list(tmp_path.glob("*.outputManifest.json"))
