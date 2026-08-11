"""Bootstrap confidence intervals for per-amplicon editing rates."""

from dataclasses import dataclass
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


@dataclass(frozen=True)
class EditingRateCIConfig:
    """Configuration for editing-rate bootstrap confidence intervals."""

    enabled: bool = False
    bootstrap_iterations: int = 10_000
    confidence_level: float = 0.95
    seed: int = 42
    batch_size: int = 64


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
