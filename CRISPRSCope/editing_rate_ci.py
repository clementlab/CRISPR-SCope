"""Bootstrap confidence intervals for per-amplicon editing rates."""

from dataclasses import dataclass
import math
import multiprocessing as mp
from typing import Dict, Iterable, List, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


RESULT_COLUMNS = [
    "amplicon",
    "all_estimate_pct",
    "all_ci_lower_pct",
    "all_ci_upper_pct",
    "all_n_cells",
    "hq_estimate_pct",
    "hq_ci_lower_pct",
    "hq_ci_upper_pct",
    "hq_n_cells",
    "hq_minus_all_pct",
    "delta_ci_lower_pct",
    "delta_ci_upper_pct",
    "valid_all_bootstrap_replicates",
    "valid_hq_bootstrap_replicates",
    "valid_delta_bootstrap_replicates",
    "bootstrap_iterations",
    "confidence_level",
    "seed",
    "status",
]


DEPTH_STABILITY_COLUMNS = [
    "amplicon",
    "cohort",
    "sample_percent",
    "sample_n_cells",
    "eligible_n_cells",
    "full_estimate_pct",
    "subsample_median_pct",
    "subsample_interval_lower_pct",
    "subsample_interval_upper_pct",
    "median_deviation_pp",
    "deviation_interval_lower_pp",
    "deviation_interval_upper_pp",
    "median_abs_deviation_from_full_pp",
    "p95_abs_deviation_from_full_pp",
    "valid_subsamples",
    "requested_iterations",
    "confidence_level",
    "seed",
    "is_full_reference",
    "status",
    "hq_full_estimate_pct",
    "relative_min_hq_edit_pct",
    "relative_plot_eligible",
    "median_relative_deviation_pct",
    "relative_deviation_interval_lower_pct",
    "relative_deviation_interval_upper_pct",
    "median_absolute_relative_deviation_pct",
    "p95_absolute_relative_deviation_pct",
]


@dataclass(frozen=True)
class EditingRateCIConfig:
    """Configuration for editing-rate bootstrap confidence intervals."""

    enabled: bool = True
    bootstrap_iterations: int = 10_000
    confidence_level: float = 0.95
    seed: int = 42
    batch_size: int = 64


@dataclass(frozen=True)
class EditingRateDepthStabilityConfig:
    """Configuration for finite-cohort editing-rate downsampling."""

    enabled: bool = True
    iterations: int = 1_000
    percentages: Tuple[float, ...] = (10.0, 25.0, 50.0, 75.0, 90.0)
    confidence_level: float = 0.95
    seed: int = 42
    relative_min_hq_edit_pct: float = 1.0


def _point_estimate(values: np.ndarray, valid: np.ndarray) -> float:
    if not np.any(valid):
        return np.nan
    return float(np.mean(values[valid]))


def _percentile_interval(values: Sequence[float], confidence_level: float) -> Tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan, np.nan
    alpha = 1.0 - confidence_level
    lower, upper = np.quantile(values, [alpha / 2.0, 1.0 - alpha / 2.0])
    return float(lower), float(upper)


def _bootstrap_means(
    values: np.ndarray,
    iterations: int,
    rng: np.random.Generator,
    batch_size: int,
) -> np.ndarray:
    """Bootstrap a mean using a fixed number of draws from eligible values."""
    values = np.asarray(values, dtype=float)
    n_values = len(values)
    bootstrap = np.empty(iterations, dtype=float)
    completed = 0
    while completed < iterations:
        this_batch_size = min(batch_size, iterations - completed)
        sampled_indices = rng.integers(0, n_values, size=(this_batch_size, n_values))
        bootstrap[completed:completed + this_batch_size] = np.take(values, sampled_indices).mean(axis=1)
        completed += this_batch_size
    return bootstrap


def _bootstrap_stratified_delta(
    hq_values: np.ndarray,
    non_hq_values: np.ndarray,
    iterations: int,
    rng: np.random.Generator,
    batch_size: int,
) -> np.ndarray:
    """Bootstrap HQ-minus-all while fixing observed HQ/non-HQ sample sizes."""
    hq_values = np.asarray(hq_values, dtype=float)
    non_hq_values = np.asarray(non_hq_values, dtype=float)
    n_hq = len(hq_values)
    n_non_hq = len(non_hq_values)
    n_all = n_hq + n_non_hq
    bootstrap = np.empty(iterations, dtype=float)
    completed = 0
    while completed < iterations:
        this_batch_size = min(batch_size, iterations - completed)
        hq_indices = rng.integers(0, n_hq, size=(this_batch_size, n_hq))
        hq_sums = np.take(hq_values, hq_indices).sum(axis=1)
        if n_non_hq:
            non_hq_indices = rng.integers(0, n_non_hq, size=(this_batch_size, n_non_hq))
            non_hq_sums = np.take(non_hq_values, non_hq_indices).sum(axis=1)
        else:
            non_hq_sums = np.zeros(this_batch_size, dtype=float)
        hq_means = hq_sums / n_hq
        all_means = (hq_sums + non_hq_sums) / n_all
        bootstrap[completed:completed + this_batch_size] = hq_means - all_means
        completed += this_batch_size
    return bootstrap


def _bootstrap_one_amplicon(job: Dict) -> Dict:
    """Compute estimates and paired bootstrap intervals for one amplicon."""
    amplicon = job["amplicon"]
    amplicon_index = job["amplicon_index"]
    mod_values = np.asarray(job["mod_values"], dtype=float)
    count_values = np.asarray(job["count_values"], dtype=float)
    hq_mask = np.asarray(job["hq_mask"], dtype=bool)
    min_reads = job["min_reads"]
    iterations = job["iterations"]
    confidence_level = job["confidence_level"]
    base_seed = job["seed"]
    batch_size = job["batch_size"]

    valid_all = np.isfinite(mod_values) & np.isfinite(count_values) & (count_values >= min_reads)
    valid_hq = valid_all & hq_mask
    n_all = int(valid_all.sum())
    n_hq = int(valid_hq.sum())
    all_estimate = _point_estimate(mod_values, valid_all)
    hq_estimate = _point_estimate(mod_values, valid_hq)
    delta_estimate = hq_estimate - all_estimate if np.isfinite(all_estimate) and np.isfinite(hq_estimate) else np.nan

    all_values = mod_values[valid_all]
    hq_values = mod_values[valid_hq]
    non_hq_values = mod_values[valid_all & ~hq_mask]
    seed_sequence = np.random.SeedSequence([base_seed, amplicon_index])
    all_seed, hq_seed, delta_seed = seed_sequence.spawn(3)

    all_bootstrap = np.array([], dtype=float)
    if n_all >= 2:
        all_bootstrap = _bootstrap_means(
            all_values, iterations, np.random.default_rng(all_seed), batch_size
        )

    hq_bootstrap = np.array([], dtype=float)
    if n_hq >= 2:
        hq_bootstrap = _bootstrap_means(
            hq_values, iterations, np.random.default_rng(hq_seed), batch_size
        )

    delta_bootstrap = np.array([], dtype=float)
    if n_all >= 2 and n_hq >= 2:
        delta_bootstrap = _bootstrap_stratified_delta(
            hq_values,
            non_hq_values,
            iterations,
            np.random.default_rng(delta_seed),
            batch_size,
        )

    all_lower, all_upper = _percentile_interval(all_bootstrap, confidence_level)
    hq_lower, hq_upper = _percentile_interval(hq_bootstrap, confidence_level)
    delta_lower, delta_upper = _percentile_interval(delta_bootstrap, confidence_level)

    status_parts = []
    if n_all < 2:
        status_parts.append("insufficient_all_cells")
    if n_hq < 2:
        status_parts.append("insufficient_hq_cells")

    return {
        "amplicon": amplicon,
        "all_estimate_pct": all_estimate,
        "all_ci_lower_pct": all_lower,
        "all_ci_upper_pct": all_upper,
        "all_n_cells": n_all,
        "hq_estimate_pct": hq_estimate,
        "hq_ci_lower_pct": hq_lower,
        "hq_ci_upper_pct": hq_upper,
        "hq_n_cells": n_hq,
        "hq_minus_all_pct": delta_estimate,
        "delta_ci_lower_pct": delta_lower,
        "delta_ci_upper_pct": delta_upper,
        "valid_all_bootstrap_replicates": len(all_bootstrap),
        "valid_hq_bootstrap_replicates": len(hq_bootstrap),
        "valid_delta_bootstrap_replicates": len(delta_bootstrap),
        "bootstrap_iterations": iterations,
        "confidence_level": confidence_level,
        "seed": base_seed,
        "status": ";".join(status_parts) if status_parts else "ok",
    }


def compute_editing_rate_confidence_intervals(
    editing_summary: pd.DataFrame,
    quality_scores: pd.DataFrame,
    high_quality_codes: Iterable[str],
    min_reads_per_amplicon_per_cell: int,
    config: EditingRateCIConfig,
    n_processes: int = 1,
) -> pd.DataFrame:
    """Compute pointwise editing-rate intervals for all and selected HQ cells."""
    if "Color" not in quality_scores.columns:
        raise ValueError("Amplicon score table must contain a 'Color' column")

    mod_columns = [column for column in editing_summary.columns if column.startswith("modPct.")]
    hq_mask = quality_scores.reindex(editing_summary.index)["Color"].isin(set(high_quality_codes)).to_numpy()
    jobs = []
    for amplicon_index, mod_column in enumerate(mod_columns):
        amplicon = mod_column.split(".", 1)[1]
        count_column = f"totCount.{amplicon}"
        if count_column not in editing_summary.columns:
            raise ValueError(f"Missing matching count column for amplicon {amplicon!r}: {count_column}")
        jobs.append(
            {
                "amplicon": amplicon,
                "amplicon_index": amplicon_index,
                "mod_values": pd.to_numeric(editing_summary[mod_column], errors="coerce").to_numpy(dtype=float),
                "count_values": pd.to_numeric(editing_summary[count_column], errors="coerce").to_numpy(dtype=float),
                "hq_mask": hq_mask,
                "min_reads": min_reads_per_amplicon_per_cell,
                "iterations": config.bootstrap_iterations,
                "confidence_level": config.confidence_level,
                "seed": config.seed,
                "batch_size": config.batch_size,
            }
        )

    if not jobs:
        return pd.DataFrame(columns=RESULT_COLUMNS)

    worker_count = min(max(1, int(n_processes)), len(jobs))
    if worker_count > 1:
        with mp.Pool(worker_count) as pool:
            results = pool.map(_bootstrap_one_amplicon, jobs)
    else:
        results = [_bootstrap_one_amplicon(job) for job in jobs]
    return pd.DataFrame(results, columns=RESULT_COLUMNS)


def _nested_subsample_means(
    values: np.ndarray,
    sample_sizes: Sequence[int],
    iterations: int,
    rng: np.random.Generator,
) -> Dict[int, np.ndarray]:
    """Calculate means from nested samples drawn without replacement."""
    values = np.asarray(values, dtype=float)
    unique_sizes = np.asarray(sorted(set(int(size) for size in sample_sizes)), dtype=int)
    if unique_sizes.size == 0:
        return {}
    if unique_sizes[0] < 1 or unique_sizes[-1] > len(values):
        raise ValueError("Subsample sizes must be between 1 and the number of values")

    sampled_means = np.empty((iterations, len(unique_sizes)), dtype=float)
    segment_starts = np.concatenate(([0], unique_sizes[:-1]))
    maximum_size = int(unique_sizes[-1])
    for iteration in range(iterations):
        sampled_indices = rng.choice(
            len(values),
            size=maximum_size,
            replace=False,
            shuffle=True,
        )
        segment_sums = np.add.reduceat(values[sampled_indices], segment_starts)
        sampled_means[iteration] = np.cumsum(segment_sums) / unique_sizes
    return {
        int(sample_size): sampled_means[:, column_index]
        for column_index, sample_size in enumerate(unique_sizes)
    }


def _empty_depth_stability_row(
    amplicon: str,
    cohort: str,
    sample_percent: float,
    sample_n_cells: int,
    eligible_n_cells: int,
    full_estimate: float,
    job: Dict,
) -> Dict:
    return {
        "amplicon": amplicon,
        "cohort": cohort,
        "sample_percent": sample_percent,
        "sample_n_cells": sample_n_cells,
        "eligible_n_cells": eligible_n_cells,
        "full_estimate_pct": full_estimate,
        "subsample_median_pct": np.nan,
        "subsample_interval_lower_pct": np.nan,
        "subsample_interval_upper_pct": np.nan,
        "median_deviation_pp": np.nan,
        "deviation_interval_lower_pp": np.nan,
        "deviation_interval_upper_pp": np.nan,
        "median_abs_deviation_from_full_pp": np.nan,
        "p95_abs_deviation_from_full_pp": np.nan,
        "valid_subsamples": 0,
        "requested_iterations": job["iterations"],
        "confidence_level": job["confidence_level"],
        "seed": job["seed"],
        "is_full_reference": sample_percent == 100.0,
        "status": "insufficient_cells",
    }


def _depth_stability_one_amplicon(job: Dict) -> List[Dict]:
    """Downsample eligible all-cell and HQ values for one amplicon."""
    amplicon = job["amplicon"]
    mod_values = np.asarray(job["mod_values"], dtype=float)
    count_values = np.asarray(job["count_values"], dtype=float)
    hq_mask = np.asarray(job["hq_mask"], dtype=bool)
    valid_all = (
        np.isfinite(mod_values)
        & np.isfinite(count_values)
        & (count_values >= job["min_reads"])
    )
    percentages = tuple(float(value) for value in job["percentages"])
    amplicon_seed = np.random.SeedSequence([job["seed"], job["amplicon_index"]])
    cohort_seeds = amplicon_seed.spawn(2)
    results: List[Dict] = []

    for cohort_index, (cohort, cohort_mask) in enumerate(
        (("all", valid_all), ("hq", valid_all & hq_mask))
    ):
        values = mod_values[cohort_mask]
        eligible_n_cells = len(values)
        full_estimate = float(np.mean(values)) if eligible_n_cells else np.nan
        percentage_sizes = {
            percentage: min(
                eligible_n_cells,
                max(1, int(math.ceil(eligible_n_cells * percentage / 100.0))),
            )
            if eligible_n_cells
            else 0
            for percentage in percentages
        }

        if eligible_n_cells < 2:
            for percentage in percentages + (100.0,):
                sample_n_cells = (
                    percentage_sizes[percentage]
                    if percentage < 100.0
                    else eligible_n_cells
                )
                results.append(
                    _empty_depth_stability_row(
                        amplicon,
                        cohort,
                        percentage,
                        sample_n_cells,
                        eligible_n_cells,
                        full_estimate,
                        job,
                    )
                )
            continue

        sample_sizes = list(percentage_sizes.values())
        sampled_means = _nested_subsample_means(
            values,
            sample_sizes,
            job["iterations"],
            np.random.default_rng(cohort_seeds[cohort_index]),
        )
        for percentage in percentages:
            sample_n_cells = percentage_sizes[percentage]
            means = sampled_means[sample_n_cells]
            lower, upper = _percentile_interval(means, job["confidence_level"])
            median = float(np.median(means))
            deviations = means - full_estimate
            deviation_lower, deviation_upper = _percentile_interval(
                deviations, job["confidence_level"]
            )
            absolute_deviations = np.abs(deviations)
            results.append(
                {
                    "amplicon": amplicon,
                    "cohort": cohort,
                    "sample_percent": percentage,
                    "sample_n_cells": sample_n_cells,
                    "eligible_n_cells": eligible_n_cells,
                    "full_estimate_pct": full_estimate,
                    "subsample_median_pct": median,
                    "subsample_interval_lower_pct": lower,
                    "subsample_interval_upper_pct": upper,
                    "median_deviation_pp": median - full_estimate,
                    "deviation_interval_lower_pp": deviation_lower,
                    "deviation_interval_upper_pp": deviation_upper,
                    "median_abs_deviation_from_full_pp": float(
                        np.median(absolute_deviations)
                    ),
                    "p95_abs_deviation_from_full_pp": float(
                        np.quantile(absolute_deviations, 0.95)
                    ),
                    "valid_subsamples": len(means),
                    "requested_iterations": job["iterations"],
                    "confidence_level": job["confidence_level"],
                    "seed": job["seed"],
                    "is_full_reference": False,
                    "status": "ok",
                }
            )

        results.append(
            {
                "amplicon": amplicon,
                "cohort": cohort,
                "sample_percent": 100.0,
                "sample_n_cells": eligible_n_cells,
                "eligible_n_cells": eligible_n_cells,
                "full_estimate_pct": full_estimate,
                "subsample_median_pct": full_estimate,
                "subsample_interval_lower_pct": full_estimate,
                "subsample_interval_upper_pct": full_estimate,
                "median_deviation_pp": 0.0,
                "deviation_interval_lower_pp": 0.0,
                "deviation_interval_upper_pp": 0.0,
                "median_abs_deviation_from_full_pp": 0.0,
                "p95_abs_deviation_from_full_pp": 0.0,
                "valid_subsamples": 1,
                "requested_iterations": job["iterations"],
                "confidence_level": job["confidence_level"],
                "seed": job["seed"],
                "is_full_reference": True,
                "status": "full_reference",
            }
        )
    return results


def _add_relative_depth_stability_metrics(
    results: pd.DataFrame,
    relative_min_hq_edit_pct: float,
) -> pd.DataFrame:
    """Add HQ-thresholded relative deviations to depth-stability results."""
    results = results.copy()
    relative_columns = [
        "median_relative_deviation_pct",
        "relative_deviation_interval_lower_pct",
        "relative_deviation_interval_upper_pct",
        "median_absolute_relative_deviation_pct",
        "p95_absolute_relative_deviation_pct",
    ]
    for column in relative_columns:
        results[column] = np.nan
    results["relative_min_hq_edit_pct"] = float(relative_min_hq_edit_pct)

    hq_reference = results.loc[
        (results["cohort"] == "hq")
        & results["is_full_reference"].astype(bool)
        & (results["status"] == "full_reference"),
        ["amplicon", "full_estimate_pct"],
    ].drop_duplicates("amplicon")
    hq_estimate_by_amplicon = hq_reference.set_index("amplicon")["full_estimate_pct"]
    results["hq_full_estimate_pct"] = results["amplicon"].map(hq_estimate_by_amplicon)
    results["relative_plot_eligible"] = (
        pd.to_numeric(results["hq_full_estimate_pct"], errors="coerce").ge(
            relative_min_hq_edit_pct
        )
    )

    full_estimate = pd.to_numeric(results["full_estimate_pct"], errors="coerce")
    relative_valid = (
        results["relative_plot_eligible"]
        & results["status"].isin(["ok", "full_reference"])
        & full_estimate.gt(0)
    )
    scale = 100.0 / full_estimate.loc[relative_valid]
    source_to_relative = {
        "median_deviation_pp": "median_relative_deviation_pct",
        "deviation_interval_lower_pp": "relative_deviation_interval_lower_pct",
        "deviation_interval_upper_pp": "relative_deviation_interval_upper_pct",
        "median_abs_deviation_from_full_pp": "median_absolute_relative_deviation_pct",
        "p95_abs_deviation_from_full_pp": "p95_absolute_relative_deviation_pct",
    }
    for source_column, relative_column in source_to_relative.items():
        source_values = pd.to_numeric(
            results.loc[relative_valid, source_column], errors="coerce"
        )
        results.loc[relative_valid, relative_column] = source_values * scale
    return results


def compute_editing_rate_depth_stability(
    editing_summary: pd.DataFrame,
    quality_scores: pd.DataFrame,
    high_quality_codes: Iterable[str],
    min_reads_per_amplicon_per_cell: int,
    config: EditingRateDepthStabilityConfig,
    n_processes: int = 1,
) -> pd.DataFrame:
    """Measure editing-rate stability under finite-cohort downsampling."""
    if "Color" not in quality_scores.columns:
        raise ValueError("Amplicon score table must contain a 'Color' column")
    relative_min_hq_edit_pct = float(config.relative_min_hq_edit_pct)
    if (
        not np.isfinite(relative_min_hq_edit_pct)
        or relative_min_hq_edit_pct <= 0
        or relative_min_hq_edit_pct > 100
    ):
        raise ValueError(
            "relative_min_hq_edit_pct must be greater than 0 and no greater than 100"
        )

    row_labels = np.asarray([str(value) for value in editing_summary.index])
    row_order = np.argsort(row_labels, kind="stable")
    ordered_index = editing_summary.index[row_order]
    ordered_summary = editing_summary.iloc[row_order]
    hq_mask = (
        quality_scores.reindex(ordered_index)["Color"]
        .isin(set(high_quality_codes))
        .to_numpy()
    )
    mod_columns = [
        column for column in ordered_summary.columns if column.startswith("modPct.")
    ]
    jobs = []
    for amplicon_index, mod_column in enumerate(mod_columns):
        amplicon = mod_column.split(".", 1)[1]
        count_column = f"totCount.{amplicon}"
        if count_column not in ordered_summary.columns:
            raise ValueError(
                f"Missing matching count column for amplicon {amplicon!r}: {count_column}"
            )
        jobs.append(
            {
                "amplicon": amplicon,
                "amplicon_index": amplicon_index,
                "mod_values": pd.to_numeric(
                    ordered_summary[mod_column], errors="coerce"
                ).to_numpy(dtype=float),
                "count_values": pd.to_numeric(
                    ordered_summary[count_column], errors="coerce"
                ).to_numpy(dtype=float),
                "hq_mask": hq_mask,
                "min_reads": min_reads_per_amplicon_per_cell,
                "iterations": config.iterations,
                "percentages": config.percentages,
                "confidence_level": config.confidence_level,
                "seed": config.seed,
            }
        )

    if not jobs:
        return pd.DataFrame(columns=DEPTH_STABILITY_COLUMNS)

    worker_count = min(max(1, int(n_processes)), len(jobs))
    if worker_count > 1:
        with mp.Pool(worker_count) as pool:
            nested_results = pool.map(_depth_stability_one_amplicon, jobs)
    else:
        nested_results = [_depth_stability_one_amplicon(job) for job in jobs]
    results = [row for amplicon_rows in nested_results for row in amplicon_rows]
    result_frame = pd.DataFrame(results)
    result_frame = _add_relative_depth_stability_metrics(
        result_frame,
        relative_min_hq_edit_pct,
    )
    return result_frame.reindex(columns=DEPTH_STABILITY_COLUMNS)


def _finite_interval_rows(results: pd.DataFrame, prefix: str) -> pd.Series:
    return (
        pd.to_numeric(results[f"{prefix}_estimate_pct"], errors="coerce").notna()
        & pd.to_numeric(results[f"{prefix}_ci_lower_pct"], errors="coerce").notna()
        & pd.to_numeric(results[f"{prefix}_ci_upper_pct"], errors="coerce").notna()
    )


def write_editing_rate_ci_plots(results: pd.DataFrame, output_root: str) -> List[Dict[str, str]]:
    """Write HQ/all and paired-delta forest plots and return report metadata."""
    plot_metadata: List[Dict[str, str]] = []
    if results.empty:
        return plot_metadata

    ordered = results.copy()
    ordered["all_estimate_pct"] = pd.to_numeric(ordered["all_estimate_pct"], errors="coerce")
    ordered = ordered.sort_values("all_estimate_pct", ascending=True, na_position="first").reset_index(drop=True)

    all_valid = _finite_interval_rows(ordered, "all")
    hq_valid = _finite_interval_rows(ordered, "hq")
    if all_valid.any() or hq_valid.any():
        fig_height = max(6.0, 0.4 * len(ordered) + 2.0)
        fig, ax = plt.subplots(figsize=(12, fig_height))
        y_positions = np.arange(len(ordered), dtype=float)
        for valid, prefix, label, color, offset in [
            (all_valid, "all", "All analyzable cells", "#4C78A8", -0.12),
            (hq_valid, "hq", "Configured high-quality cells", "#F58518", 0.12),
        ]:
            subset = ordered.loc[valid]
            positions = y_positions[valid.to_numpy()] + offset
            estimate = subset[f"{prefix}_estimate_pct"].to_numpy(dtype=float)
            lower = subset[f"{prefix}_ci_lower_pct"].to_numpy(dtype=float)
            upper = subset[f"{prefix}_ci_upper_pct"].to_numpy(dtype=float)
            ax.errorbar(
                estimate,
                positions,
                xerr=np.vstack([estimate - lower, upper - estimate]),
                fmt="o",
                capsize=3,
                label=label,
                color=color,
            )
        ax.set_yticks(y_positions)
        ax.set_yticklabels(ordered["amplicon"])
        ax.set_xlim(0, 100)
        ax.set_xlabel("Mean inferred allele editing percentage")
        ax.set_ylabel("Amplicon")
        ax.set_title("Amplicon Editing-Rate Confidence Intervals")
        ax.legend()
        ax.grid(axis="x", alpha=0.25)
        fig.tight_layout()
        plot_root = output_root + ".10_EditingRateConfidenceIntervals"
        fig.savefig(plot_root + ".pdf", bbox_inches="tight")
        fig.savefig(plot_root + ".png", bbox_inches="tight")
        plt.close(fig)
        plot_metadata.append(
            {
                "plot_name": plot_root,
                "plot_title": "Amplicon editing-rate confidence intervals",
                "plot_label": (
                    "Pointwise bootstrap confidence intervals for all analyzable "
                    "and configured high-quality cells."
                ),
            }
        )

    delta_valid = (
        pd.to_numeric(ordered["hq_minus_all_pct"], errors="coerce").notna()
        & pd.to_numeric(ordered["delta_ci_lower_pct"], errors="coerce").notna()
        & pd.to_numeric(ordered["delta_ci_upper_pct"], errors="coerce").notna()
    )
    if delta_valid.any():
        fig_height = max(6.0, 0.4 * len(ordered) + 2.0)
        fig, ax = plt.subplots(figsize=(12, fig_height))
        subset = ordered.loc[delta_valid]
        positions = np.arange(len(ordered), dtype=float)[delta_valid.to_numpy()]
        estimate = subset["hq_minus_all_pct"].to_numpy(dtype=float)
        lower = subset["delta_ci_lower_pct"].to_numpy(dtype=float)
        upper = subset["delta_ci_upper_pct"].to_numpy(dtype=float)
        ax.errorbar(
            estimate,
            positions,
            xerr=np.vstack([estimate - lower, upper - estimate]),
            fmt="o",
            capsize=3,
            color="#54A24B",
        )
        ax.axvline(0, color="black", linestyle="--", linewidth=1)
        ax.set_yticks(np.arange(len(ordered), dtype=float))
        ax.set_yticklabels(ordered["amplicon"])
        ax.set_xlim(-100, 100)
        ax.set_xlabel("HQ minus all-cell editing percentage points")
        ax.set_ylabel("Amplicon")
        ax.set_title("Editing-Rate Sensitivity to Cell-Quality Selection")
        ax.grid(axis="x", alpha=0.25)
        fig.tight_layout()
        plot_root = output_root + ".11_EditingRateQualityDelta"
        fig.savefig(plot_root + ".pdf", bbox_inches="tight")
        fig.savefig(plot_root + ".png", bbox_inches="tight")
        plt.close(fig)
        plot_metadata.append(
            {
                "plot_name": plot_root,
                "plot_title": "Editing-rate quality-selection sensitivity",
                "plot_label": (
                    "Paired bootstrap interval for the configured high-quality "
                    "estimate minus the all-analyzable-cell estimate."
                ),
            }
        )
    return plot_metadata


def _depth_stability_amplicon_order(results: pd.DataFrame) -> List[str]:
    """Order amplicons by HQ, then all-cell, full-cohort editing estimates."""
    references = results.loc[
        results["is_full_reference"].astype(bool)
        & (results["status"] == "full_reference"),
        ["amplicon", "cohort", "full_estimate_pct"],
    ].copy()
    references["full_estimate_pct"] = pd.to_numeric(
        references["full_estimate_pct"], errors="coerce"
    )
    reference_lookup = references.pivot_table(
        index="amplicon",
        columns="cohort",
        values="full_estimate_pct",
        aggfunc="first",
    )

    def sort_key(amplicon: str) -> Tuple:
        if amplicon in reference_lookup.index:
            hq_estimate = (
                reference_lookup.at[amplicon, "hq"]
                if "hq" in reference_lookup
                else np.nan
            )
            all_estimate = (
                reference_lookup.at[amplicon, "all"]
                if "all" in reference_lookup
                else np.nan
            )
        else:
            hq_estimate = np.nan
            all_estimate = np.nan
        hq_is_missing = not np.isfinite(hq_estimate)
        all_is_missing = not np.isfinite(all_estimate)
        return (
            hq_is_missing,
            -float(hq_estimate) if not hq_is_missing else 0.0,
            all_is_missing,
            -float(all_estimate) if not all_is_missing else 0.0,
            str(amplicon),
        )

    amplicons = results["amplicon"].drop_duplicates().tolist()
    return sorted(amplicons, key=sort_key)


def _create_depth_stability_figure(
    plotted: pd.DataFrame,
    amplicon_order: Sequence[str],
    median_column: str,
    lower_column: str,
    upper_column: str,
    title: str,
    y_label: str,
    subtitle: str = "",
):
    """Create a consistently labeled three-column stability facet grid."""
    n_columns = 3
    n_rows = int(math.ceil(len(amplicon_order) / n_columns))
    fig, axes = plt.subplots(
        n_rows,
        n_columns,
        figsize=(18, max(5.0, 3.8 * n_rows)),
        sharex=True,
        sharey=True,
        squeeze=False,
    )
    axes_flat = axes.ravel()
    finite_bounds = plotted[[lower_column, upper_column]].to_numpy(dtype=float)
    maximum_absolute_bound = float(np.nanmax(np.abs(finite_bounds)))
    y_limit = max(0.1, maximum_absolute_bound * 1.1)
    styles = {
        "all": ("All analyzable cells", "#4C78A8"),
        "hq": ("Configured high-quality cells", "#F58518"),
    }
    legend_handles = []
    legend_labels = []

    for facet_index, (axis, amplicon) in enumerate(zip(axes_flat, amplicon_order)):
        amplicon_rows = plotted.loc[plotted["amplicon"] == amplicon]
        for cohort, (label, color) in styles.items():
            cohort_rows = amplicon_rows.loc[
                amplicon_rows["cohort"] == cohort
            ].sort_values("sample_percent")
            if cohort_rows.empty:
                continue
            percentages = cohort_rows["sample_percent"].to_numpy(dtype=float)
            medians = cohort_rows[median_column].to_numpy(dtype=float)
            lower = cohort_rows[lower_column].to_numpy(dtype=float)
            upper = cohort_rows[upper_column].to_numpy(dtype=float)
            line = axis.plot(
                percentages,
                medians,
                marker="o",
                markersize=3,
                linewidth=1.3,
                color=color,
                label=label,
            )[0]
            axis.fill_between(percentages, lower, upper, color=color, alpha=0.2)
            if label not in legend_labels:
                legend_handles.append(line)
                legend_labels.append(label)
        axis.axhline(0, color="black", linestyle="--", linewidth=0.8, alpha=0.7)
        axis.set_title(amplicon, fontsize=8)
        axis.set_xlim(0, 102)
        axis.set_ylim(-y_limit, y_limit)
        axis.grid(alpha=0.2)
        axis.tick_params(
            axis="x",
            labelbottom=True,
            labelsize=7,
        )
        axis.tick_params(
            axis="y",
            labelleft=(facet_index % n_columns == 0),
            labelsize=7,
        )

    for axis in axes_flat[len(amplicon_order):]:
        axis.set_visible(False)

    tick_values = sorted(plotted["sample_percent"].unique())
    for axis in axes_flat[:len(amplicon_order)]:
        axis.set_xticks(tick_values)

    x_label = "Eligible cohort retained (%)"
    for row_index in range(n_rows):
        row_start = row_index * n_columns
        populated_count = min(n_columns, len(amplicon_order) - row_start)
        if populated_count <= 0:
            continue
        label_column = 1 if populated_count >= 2 else 0
        axes[row_index, label_column].set_xlabel(x_label, fontsize=9, labelpad=5)
        axes[row_index, 0].set_ylabel(y_label, fontsize=8, labelpad=7)

    fig.suptitle(title, fontsize=16, y=0.997)
    if subtitle:
        fig.text(0.5, 0.981, subtitle, ha="center", va="top", fontsize=10)
        legend_y = 0.967
        layout_top = 0.94
    else:
        legend_y = 0.981
        layout_top = 0.955
    if legend_handles:
        fig.legend(
            legend_handles,
            legend_labels,
            loc="upper center",
            bbox_to_anchor=(0.5, legend_y),
            ncol=2,
        )
    fig.tight_layout(rect=(0.02, 0.01, 1.0, layout_top), h_pad=2.4, w_pad=1.0)
    return fig


def _save_depth_stability_figure(fig, plot_root: str) -> None:
    fig.savefig(plot_root + ".pdf", bbox_inches="tight")
    fig.savefig(plot_root + ".png", bbox_inches="tight")
    plt.close(fig)


def write_editing_rate_depth_stability_plot(
    results: pd.DataFrame,
    output_root: str,
) -> List[Dict[str, str]]:
    """Write absolute and thresholded-relative editing-rate stability plots."""
    if results.empty:
        return []

    plotted = results.loc[
        results["status"].isin(["ok", "full_reference"])
    ].copy()
    numeric_columns = [
        "sample_percent",
        "full_estimate_pct",
        "median_deviation_pp",
        "deviation_interval_lower_pp",
        "deviation_interval_upper_pp",
        "median_relative_deviation_pct",
        "relative_deviation_interval_lower_pct",
        "relative_deviation_interval_upper_pct",
    ]
    for column in numeric_columns:
        plotted[column] = pd.to_numeric(plotted[column], errors="coerce")
    absolute_plotted = plotted.dropna(
        subset=[
            "sample_percent",
            "median_deviation_pp",
            "deviation_interval_lower_pp",
            "deviation_interval_upper_pp",
        ]
    )
    if absolute_plotted.empty:
        return []

    amplicon_order = _depth_stability_amplicon_order(absolute_plotted)
    absolute_figure = _create_depth_stability_figure(
        absolute_plotted,
        amplicon_order,
        median_column="median_deviation_pp",
        lower_column="deviation_interval_lower_pp",
        upper_column="deviation_interval_upper_pp",
        title="Editing-Rate Stability Across Cell-Depth Downsampling",
        y_label="Deviation from full-cohort editing rate (percentage points)",
    )
    absolute_plot_root = output_root + ".12_EditingRateDepthStability"
    _save_depth_stability_figure(absolute_figure, absolute_plot_root)
    plot_metadata = [
        {
            "plot_name": absolute_plot_root,
            "plot_title": "Editing-rate cell-depth stability",
            "plot_label": (
                "Finite-cohort downsampling stability bands for all analyzable "
                "and configured high-quality cells. Bands show sensitivity to "
                "retained cell depth within this run, not uncertainty across "
                "biological replicates."
            ),
        }
    ]

    relative_eligible = plotted["relative_plot_eligible"].fillna(False).astype(bool)
    relative_plotted = plotted.loc[relative_eligible].dropna(
        subset=[
            "sample_percent",
            "median_relative_deviation_pct",
            "relative_deviation_interval_lower_pct",
            "relative_deviation_interval_upper_pct",
        ]
    )
    if relative_plotted.empty:
        return plot_metadata

    relative_amplicons = set(relative_plotted["amplicon"])
    relative_order = [
        amplicon for amplicon in amplicon_order if amplicon in relative_amplicons
    ]
    threshold = float(
        pd.to_numeric(
            relative_plotted["relative_min_hq_edit_pct"], errors="coerce"
        ).dropna().iloc[0]
    )
    relative_figure = _create_depth_stability_figure(
        relative_plotted,
        relative_order,
        median_column="median_relative_deviation_pct",
        lower_column="relative_deviation_interval_lower_pct",
        upper_column="relative_deviation_interval_upper_pct",
        title="Relative Editing-Rate Stability Across Cell-Depth Downsampling",
        subtitle=(
            "Amplicons included when the full configured-HQ editing rate is "
            f"at least {threshold:g}%"
        ),
        y_label="Relative deviation from full-cohort editing rate (%)",
    )
    relative_plot_root = output_root + ".13_EditingRateRelativeDepthStability"
    _save_depth_stability_figure(relative_figure, relative_plot_root)
    plot_metadata.append(
        {
            "plot_name": relative_plot_root,
            "plot_title": "Relative editing-rate cell-depth stability",
            "plot_label": (
                "Relative finite-cohort downsampling stability bands for "
                "amplicons meeting the configured full-HQ editing-rate threshold. "
                "Each cohort is normalized to its own full-cohort editing rate."
            ),
        }
    )
    return plot_metadata
