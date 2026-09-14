"""Declarative names, metadata, and lifecycle tracking for run artifacts.

This module deliberately preserves CRISPRSCope's established output names.  It
centralizes the description of reportable artifacts without changing the
existing stage-file naming API used by FASTQ/BAM intermediates.
"""

from __future__ import annotations

import json
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Mapping, Sequence


@dataclass(frozen=True)
class OutputSpec:
    """One fixed, run-level output artifact."""

    key: str
    suffix: str
    kind: str
    formats: tuple[str, ...] = ()
    title: str | None = None
    label: str | None = None
    data_links: tuple[tuple[str, str], ...] = ()
    optional: bool = False


@dataclass(frozen=True)
class ArtifactFamily:
    """A parameterized intermediate-artifact family.

    ``resolve`` accepts the existing ``build_stage_filename`` callable.  This
    keeps intermediate naming in its established implementation while making
    the family available to future cache/output work.
    """

    key: str
    directory_suffix: str
    description: str

    def directory(self, output_root: str) -> str:
        return output_root + self.directory_suffix

    def resolve(self, filename_builder: Callable[..., str], **kwargs: object) -> str:
        return filename_builder(**kwargs)


ARTIFACT_SPECS: tuple[OutputSpec, ...] = (
    OutputSpec("editing_summary", ".editingSummary.txt", "table"),
    OutputSpec("filtered_editing_summary", ".filteredEditingSummary.txt", "table"),
    OutputSpec("editing_summary_pseudobulk", ".editingSummaryPseudobulk.txt", "table"),
    OutputSpec("filtered_editing_summary_pseudobulk", ".filteredEditingSummaryPseudobulk.txt", "table"),
    OutputSpec("amplicon_score", ".amplicon_score.txt", "table"),
    OutputSpec("valid_amplicons", ".splitReads.valid_amps.txt", "table"),
    OutputSpec("aligned_read_counts", ".splitReads.aligned.txt", "table"),
    OutputSpec("unaligned_read_counts", ".splitReads.unaligned.txt", "table"),
    OutputSpec("amplicon_classification", ".splitReads.amp_classification.txt", "table"),
    OutputSpec("report", ".html", "report"),
    OutputSpec("h5ad", ".h5ad", "h5ad", optional=True),
    OutputSpec("editing_rate_ci", ".editingRateConfidenceIntervals.txt", "table", optional=True),
    OutputSpec("editing_rate_unconditional_permutation", ".editingRateUnconditionalPermutation.txt", "table", optional=True),
    OutputSpec("editing_rate_unconditional_simulations", ".editingRateUnconditionalPermutationSimulations.txt", "table", optional=True),
    OutputSpec("editing_rate_depth_stability", ".editingRateDepthStability.txt", "table", optional=True),
    OutputSpec("output_manifest", ".outputManifest.json", "manifest", optional=True),
    OutputSpec(
        "log_log_plot", ".01_Log-Log", "plot", ("png", "pdf"),
        "Log-Log Plot",
        "Log scale plot of the Barcode Rank (X) vs. Read Count (Y) with color coding according to Amplicon Score category.",
    ),
    OutputSpec(
        "log_log_filtered_plot", ".01_Log-Log_filtered", "plot", ("png", "pdf"),
        "Log-Log Plot",
        "Log scale plot of the Barcode Rank (X) vs. Read Count (Y) with color coding according to Amplicon Score category.",
    ),
    OutputSpec(
        "cell_count_per_amplicon_plot", ".02_CellCountPerAmplicon_filtered", "plot", ("png", "pdf"),
        "Cell count per amplicon with minimum specified coverage",
        "Plotting the number of cells covering an amplicon at a given read cutoff. The log cell counts are on the Y-axis and the amplicon is on the X-axis.",
    ),
    OutputSpec(
        "amplicon_covered_per_cell_plot", ".03_AmpliconCoveredPerCell_filtered", "plot", ("png", "pdf"),
        "Amplicons covered per cell with minimum specified coverage",
        "The number of cells with amplicon coverage at a given read cutoff. The number of amplicons covered is on the X-axis and the cell count is on the Y-axis.",
    ),
    OutputSpec(
        "modification_percentage_plot", ".04_ModPercentagePerAmp_filtered", "plot", ("png", "pdf"),
        "Average modification by target",
        "The average modification percentage (Y) of an amplicon (X) with a minimum specified read coverage.",
    ),
    OutputSpec(
        "amplicon_score_plot", ".05_Amplicon_Score", "plot", ("png", "pdf"),
        "Supported Amplicon Breadth Plot",
        "Barcode rank versus the fraction of usable amplicons meeting the configured read-support threshold. Dashed lines show the configured breadth and barcode-rank classification boundaries.",
        (("Amplicon Score", "amplicon_score"),),
    ),
    OutputSpec(
        "edit_combinations_plot", ".06_EditCombinations", "plot", ("png", "pdf"),
        "Editing Sites and Intersections",
        "An upset plot that displays the most common edits and edit combinations.",
        (("Filtered modification percentages (modPct)", "filtered_editing_summary"), ("Filtered cell barcodes", "amplicon_score")),
    ),
    OutputSpec(
        "edit_histogram_plot", ".07_EditHistogram", "plot", ("png", "pdf"),
        "Editing Count Histogram",
        "A histogram that displays the number of edited sites in each barcode.",
        (("Filtered modification percentages (modPct)", "filtered_editing_summary"), ("Filtered cell barcodes", "amplicon_score")),
    ),
    OutputSpec(
        "cell_coverage_plot", ".08_CellCoverage", "plot", ("png", "pdf"),
        "Read Counts per Barcode",
        "A barplot displaying the difference in average read count per barcode in high quality and low quality cells.",
        (("Barcode Read Counts (totCounts)", "editing_summary_pseudobulk"), ("Filtered cell barcodes", "amplicon_score")),
    ),
    OutputSpec(
        "amplicon_coverage_plot", ".09_AmpliconCoverage", "plot", ("png", "pdf"),
        "Average Amplicon Read Coverage",
        "The average read counts covering each amplicon.",
        (("Amplicon Read Counts (totCounts)", "editing_summary_pseudobulk"), ("Filtered cell barcodes", "amplicon_score")),
    ),
    OutputSpec(
        "cell_coverage_boxplot", ".09_CellCoverageBoxplot", "plot", ("png", "pdf"),
        "Read Counts per Barcode Boxplot",
        "A boxplot displaying the read counts per barcode in high quality and low quality cells.",
        (("Barcode Read Counts (totCounts)", "editing_summary_pseudobulk"), ("Filtered cell barcodes", "amplicon_score")),
    ),
    OutputSpec(
        "editing_rate_confidence_intervals_plot", ".10_EditingRateConfidenceIntervals", "plot", ("png", "pdf"),
        "Amplicon editing-rate confidence intervals",
        "Pointwise bootstrap confidence intervals for all analyzable and configured analysis-group cells among amplicons with a coverage-adjusted BH p-value at or below 0.05.",
        (("Editing-rate confidence intervals", "editing_rate_ci"),), True,
    ),
    OutputSpec(
        "editing_rate_coverage_adjusted_effects_plot", ".11_EditingRateCoverageAdjustedEffects", "plot", ("png", "pdf"),
        "Raw vs. read-depth adjusted editing rate",
        "Observed InGroup-versus-OutGroup editing-rate differences before and after read-depth adjustment.",
        (("Editing-rate confidence intervals", "editing_rate_ci"),), True,
    ),
    OutputSpec(
        "editing_rate_unconditional_permutation_plot", ".12_EditingRateUnconditionalPermutation", "plot", ("png", "pdf"),
        "Unconditional editing-rate permutation distribution",
        "Permutation distribution of InGroup-sized subset means relative to the observed InGroup mean.",
        (("Unconditional permutation summary", "editing_rate_unconditional_permutation"), ("Unconditional permutation simulations", "editing_rate_unconditional_simulations")), True,
    ),
    OutputSpec(
        "editing_rate_observed_centered_permutation_swarm_plot", ".13_EditingRateObservedCenteredPermutationSwarm", "plot", ("png", "pdf"),
        "Observed-centered unconditional permutation swarm",
        "Simulated InGroup-sized subset means relative to the observed InGroup editing rate for each amplicon; zero marks the observed rate.",
        (("Unconditional permutation summary", "editing_rate_unconditional_permutation"), ("Unconditional permutation simulations", "editing_rate_unconditional_simulations")), True,
    ),
    OutputSpec(
        "editing_rate_depth_stability_plot", ".14_EditingRateDepthStability", "plot", ("png", "pdf"),
        "Editing-rate cell-depth stability",
        "Finite-cohort downsampling stability bands for all analyzable and configured analysis-group cells. Bands show sensitivity to retained cell depth within this run, not uncertainty across biological replicates.",
        (("Editing-rate depth stability", "editing_rate_depth_stability"),), True,
    ),
)

INTERMEDIATE_FAMILIES: tuple[ArtifactFamily, ...] = (
    ArtifactFamily("amplicon_fastq", ".seq_by_amplicon", "Per-amplicon staged FASTQ files"),
    ArtifactFamily("crispresso", ".crispresso", "Per-amplicon CRISPResso output directories"),
    ArtifactFamily("filtered_crispresso", ".crispresso.filtered", "Filtered per-amplicon CRISPResso output directories"),
)


class OutputContext:
    """Resolve registered paths while preserving established naming behavior."""

    def __init__(
        self,
        output_root: str,
        h5ad_output: str | None = None,
        specs: Sequence[OutputSpec] = ARTIFACT_SPECS,
    ):
        self.output_root = str(output_root)
        self.h5ad_output = str(h5ad_output) if h5ad_output else None
        self._specs = {spec.key: spec for spec in specs}
        if len(self._specs) != len(specs):
            raise RuntimeError("Output artifact keys must be unique")

    def spec(self, key: str) -> OutputSpec:
        try:
            return self._specs[key]
        except KeyError as error:
            raise KeyError(f"Unknown output artifact key: {key}") from error

    def path(self, key: str) -> str:
        spec = self.spec(key)
        if key == "h5ad" and self.h5ad_output:
            return self.h5ad_output
        return self.output_root + spec.suffix

    def plot_root(self, key: str) -> str:
        spec = self.spec(key)
        if spec.kind != "plot":
            raise ValueError(f"Artifact {key!r} is not a plot")
        return self.path(key)

    def paths(self, key: str) -> tuple[str, ...]:
        spec = self.spec(key)
        base = self.path(key)
        if spec.formats:
            return tuple(f"{base}.{extension}" for extension in spec.formats)
        return (base,)

    def data_links(self, key: str) -> list[tuple[str, str]]:
        return [(label, self.path(data_key)) for label, data_key in self.spec(key).data_links]

    def plot_metadata(self, key: str) -> dict[str, object]:
        spec = self.spec(key)
        if spec.kind != "plot":
            raise ValueError(f"Artifact {key!r} is not a plot")
        return {
            "plot_name": self.plot_root(key),
            "plot_title": spec.title,
            "plot_label": spec.label,
            "plot_datas": self.data_links(key),
        }

    def family(self, key: str) -> ArtifactFamily:
        for family in INTERMEDIATE_FAMILIES:
            if family.key == key:
                return family
        raise KeyError(f"Unknown output artifact family: {key}")

    def key_for_suffix(self, suffix: str) -> str | None:
        """Return the registered key for a legacy suffix, if any."""
        for spec in ARTIFACT_SPECS:
            if spec.suffix == suffix:
                return spec.key
        return None

    def remove(self, keys: Iterable[str]) -> list[str]:
        removed: list[str] = []
        for key in keys:
            for path in self.paths(key):
                try:
                    os.remove(path)
                except FileNotFoundError:
                    continue
                removed.append(path)
        return removed


class OutputManifest:
    """Ordered lifecycle record for one pipeline run's registered artifacts."""

    SCHEMA_VERSION = 1

    def __init__(self, context: OutputContext):
        self.context = context
        self.status = "running"
        self.active_stage: str | None = None
        self.failure: dict[str, str] | None = None
        self._events: dict[str, dict[str, object]] = {
            spec.key: {"status": "not_reached"} for spec in ARTIFACT_SPECS
        }

    def set_stage(self, stage: str) -> None:
        self.active_stage = stage

    def mark_written(self, key: str) -> None:
        self.context.spec(key)
        self._events[key] = {"status": "written"}

    def mark_skipped(self, key: str, reason: str) -> None:
        self.context.spec(key)
        self._events[key] = {"status": "skipped", "reason": reason}

    def mark_removed_stale(self, key: str, reason: str = "stale output from an earlier run") -> None:
        self.context.spec(key)
        self._events[key] = {"status": "removed_stale", "reason": reason}

    def mark_existing(self) -> None:
        for spec in ARTIFACT_SPECS:
            if any(os.path.exists(path) for path in self.context.paths(spec.key)):
                if self._events[spec.key]["status"] == "not_reached":
                    self.mark_written(spec.key)

    def complete(self) -> None:
        self.mark_existing()
        self.status = "completed"
        self.active_stage = None

    def fail(self, stage: str, error: BaseException) -> None:
        self.mark_existing()
        self.status = "failed"
        self.active_stage = stage
        self.failure = {"type": type(error).__name__, "message": str(error)}

    def as_dict(self) -> dict[str, object]:
        artifacts = []
        for spec in ARTIFACT_SPECS:
            event = self._events[spec.key]
            entry: dict[str, object] = {
                "key": spec.key,
                "kind": spec.kind,
                "paths": list(self.context.paths(spec.key)),
                "status": event["status"],
            }
            if spec.data_links:
                entry["data_artifact_keys"] = [key for _, key in spec.data_links]
            if "reason" in event:
                entry["reason"] = event["reason"]
            artifacts.append(entry)
        families = []
        for family in INTERMEDIATE_FAMILIES:
            directory = family.directory(self.context.output_root)
            exists = os.path.isdir(directory)
            file_count = (
                sum(len(files) for _, _, files in os.walk(directory)) if exists else 0
            )
            families.append(
                {
                    "key": family.key,
                    "directory": directory,
                    "description": family.description,
                    "exists": exists,
                    "file_count": file_count,
                }
            )
        result: dict[str, object] = {
            "schema_version": self.SCHEMA_VERSION,
            "output_root": self.context.output_root,
            "status": self.status,
            "artifacts": artifacts,
            "artifact_families": families,
        }
        if self.failure is not None:
            result["failure"] = {"stage": self.active_stage, **self.failure}
        return result

    def write(self) -> str:
        path = Path(self.context.path("output_manifest"))
        path.parent.mkdir(parents=True, exist_ok=True)
        previous_event = self._events["output_manifest"]
        self.mark_written("output_manifest")
        payload = json.dumps(self.as_dict(), indent=2, sort_keys=False) + "\n"
        fd, temporary_path = tempfile.mkstemp(
            prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
        )
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                handle.write(payload)
            os.replace(temporary_path, path)
        except BaseException:
            self._events["output_manifest"] = previous_event
            try:
                os.unlink(temporary_path)
            except FileNotFoundError:
                pass
            raise
        return str(path)
