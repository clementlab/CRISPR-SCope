"""Bootstrap, permutation, and coverage-adjusted editing-rate analyses."""

import multiprocessing as mp
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd

from .editing_rate_common import (
    EditingRateCIConfig,
    RESULT_COLUMNS,
    UNCONDITIONAL_PERMUTATION_SIMULATION_COLUMNS,
    UNCONDITIONAL_PERMUTATION_SUMMARY_COLUMNS,
    _percentile_interval,
    _point_estimate,
    categorical_mod_pct_to_cell_edit_pct,
)

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


def _permutation_hq_minus_all(
    all_values: np.ndarray,
    n_hq: int,
    iterations: int,
    rng: np.random.Generator,
    batch_size: int,
) -> np.ndarray:
    """Permute fixed-size HQ labels and return HQ-minus-all null effects."""
    all_values = np.asarray(all_values, dtype=float)
    n_all = len(all_values)
    all_mean = all_values.mean()
    null_effects = np.empty(iterations, dtype=float)
    completed = 0
    while completed < iterations:
        this_batch_size = min(batch_size, iterations - completed)
        # Selecting the smallest random keys gives an independent, uniformly
        # sampled subset without replacement for every permutation row.
        random_keys = rng.random((this_batch_size, n_all))
        hq_indices = np.argpartition(random_keys, n_hq - 1, axis=1)[:, :n_hq]
        sampled_hq = np.take(all_values, hq_indices).mean(axis=1)
        null_effects[completed:completed + this_batch_size] = sampled_hq - all_mean
        completed += this_batch_size
    return null_effects


def _two_sided_permutation_p_value(null_effects: np.ndarray, observed_effect: float) -> float:
    """Return a finite-sample-corrected two-sided Monte Carlo p-value."""
    null_effects = np.asarray(null_effects, dtype=float)
    finite_null = null_effects[np.isfinite(null_effects)]
    if finite_null.size == 0 or not np.isfinite(observed_effect):
        return np.nan
    extreme_count = np.count_nonzero(np.abs(finite_null) >= abs(observed_effect))
    return float((extreme_count + 1) / (finite_null.size + 1))


def _benjamini_hochberg(p_values: Sequence[float]) -> np.ndarray:
    """Adjust finite p-values with the Benjamini-Hochberg FDR procedure."""
    p_values = np.asarray(p_values, dtype=float)
    adjusted = np.full(p_values.shape, np.nan, dtype=float)
    finite_indices = np.flatnonzero(np.isfinite(p_values))
    if finite_indices.size == 0:
        return adjusted

    finite_p_values = p_values[finite_indices]
    order = np.argsort(finite_p_values, kind="mergesort")
    ordered_p_values = finite_p_values[order]
    ranks = np.arange(1, len(ordered_p_values) + 1, dtype=float)
    ordered_adjusted = ordered_p_values * len(ordered_p_values) / ranks
    ordered_adjusted = np.minimum.accumulate(ordered_adjusted[::-1])[::-1]
    ordered_adjusted = np.clip(ordered_adjusted, 0.0, 1.0)

    adjusted_indices = finite_indices[order]
    adjusted[adjusted_indices] = ordered_adjusted
    return adjusted


def _coverage_bin_ids(
    count_values: np.ndarray,
    exact_max_reads: int,
    bin_width_reads: int,
) -> np.ndarray:
    """Assign exact low-depth and fixed-width higher-depth coverage bins."""
    count_values = np.asarray(count_values, dtype=float)
    rounded_counts = np.rint(count_values)
    if not np.allclose(count_values, rounded_counts, rtol=0.0, atol=1e-9):
        raise ValueError("Per-amplicon read counts must be integer-valued")
    integer_counts = rounded_counts.astype(np.int64)
    if np.any(integer_counts < 0):
        raise ValueError("Per-amplicon read counts must be non-negative")

    bin_ids = integer_counts.copy()
    binned = integer_counts > exact_max_reads
    bin_ids[binned] = (
        exact_max_reads
        + 1
        + (integer_counts[binned] - exact_max_reads - 1) // bin_width_reads
    )
    return bin_ids


def _sample_empirical_means(
    values: np.ndarray,
    sample_size: int,
    iterations: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Bootstrap means from an empirical discrete distribution."""
    unique_values, counts = np.unique(np.asarray(values, dtype=float), return_counts=True)
    probabilities = counts / counts.sum()
    sampled_counts = rng.multinomial(sample_size, probabilities, size=iterations)
    return (sampled_counts @ unique_values) / sample_size


def _sample_without_replacement_sums(
    values: np.ndarray,
    sample_size: int,
    iterations: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Sample sums without replacement from a discrete finite population."""
    unique_values, counts = np.unique(np.asarray(values, dtype=float), return_counts=True)
    sampled_counts = rng.multivariate_hypergeometric(
        counts,
        sample_size,
        size=iterations,
    )
    return sampled_counts @ unique_values


def _coverage_adjusted_resampling(
    mod_values: np.ndarray,
    count_values: np.ndarray,
    group_mask: np.ndarray,
    bootstrap_iterations: int,
    permutation_iterations: int,
    confidence_level: float,
    exact_max_reads: int,
    bin_width_reads: int,
    bootstrap_rng: np.random.Generator,
    permutation_rng: np.random.Generator,
) -> Dict:
    """Estimate and resample a group effect within comparable coverage bins."""
    mod_values = np.asarray(mod_values, dtype=float)
    count_values = np.asarray(count_values, dtype=float)
    group_mask = np.asarray(group_mask, dtype=bool)
    bin_ids = _coverage_bin_ids(
        count_values,
        exact_max_reads=exact_max_reads,
        bin_width_reads=bin_width_reads,
    )

    n_group_total = int(group_mask.sum())
    n_non_group_total = int((~group_mask).sum())
    n_group_common = 0
    n_non_group_common = 0
    mixed_bin_count = 0
    total_weight = 0.0
    observed_numerator = 0.0
    standardized_group_numerator = 0.0
    standardized_non_group_numerator = 0.0
    coverage_difference_numerator = 0.0
    bootstrap_numerator = np.zeros(bootstrap_iterations, dtype=float)
    permutation_numerator = np.zeros(permutation_iterations, dtype=float)

    for bin_id in np.unique(bin_ids):
        in_bin = bin_ids == bin_id
        bin_group = in_bin & group_mask
        bin_non_group = in_bin & ~group_mask
        n_group = int(bin_group.sum())
        n_non_group = int(bin_non_group.sum())
        if n_group == 0 or n_non_group == 0:
            continue

        mixed_bin_count += 1
        n_group_common += n_group
        n_non_group_common += n_non_group
        group_values = mod_values[bin_group]
        non_group_values = mod_values[bin_non_group]
        pooled_values = mod_values[in_bin]
        weight = n_group * n_non_group / float(n_group + n_non_group)
        total_weight += weight
        group_mean = group_values.mean()
        non_group_mean = non_group_values.mean()
        standardized_group_numerator += weight * group_mean
        standardized_non_group_numerator += weight * non_group_mean
        observed_numerator += weight * (group_mean - non_group_mean)
        coverage_difference_numerator += weight * (
            count_values[bin_group].mean() - count_values[bin_non_group].mean()
        )

        bootstrap_group_means = _sample_empirical_means(
            group_values,
            n_group,
            bootstrap_iterations,
            bootstrap_rng,
        )
        bootstrap_non_group_means = _sample_empirical_means(
            non_group_values,
            n_non_group,
            bootstrap_iterations,
            bootstrap_rng,
        )
        bootstrap_numerator += weight * (
            bootstrap_group_means - bootstrap_non_group_means
        )

        sampled_group_sums = _sample_without_replacement_sums(
            pooled_values,
            n_group,
            permutation_iterations,
            permutation_rng,
        )
        permutation_numerator += (
            sampled_group_sums - n_group * pooled_values.mean()
        )

    group_retained_pct = (
        100.0 * n_group_common / n_group_total if n_group_total else np.nan
    )
    non_group_retained_pct = (
        100.0 * n_non_group_common / n_non_group_total
        if n_non_group_total
        else np.nan
    )
    status_parts = []
    if mixed_bin_count == 0:
        status_parts.append("no_mixed_coverage_bins")
    if n_group_common < 2:
        status_parts.append("insufficient_coverage_adjusted_in_group_cells")
    if n_non_group_common < 2:
        status_parts.append("insufficient_coverage_adjusted_out_group_cells")

    effect = np.nan
    standardized_group_mean = np.nan
    standardized_non_group_mean = np.nan
    ci_lower = np.nan
    ci_upper = np.nan
    permutation_p_value = np.nan
    within_bin_coverage_difference = np.nan
    valid_bootstrap_replicates = 0
    valid_permutation_replicates = 0
    if total_weight > 0:
        standardized_group_mean = float(standardized_group_numerator / total_weight)
        standardized_non_group_mean = float(
            standardized_non_group_numerator / total_weight
        )
        effect = float(observed_numerator / total_weight)
        within_bin_coverage_difference = float(
            coverage_difference_numerator / total_weight
        )
    if total_weight > 0 and n_group_common >= 2 and n_non_group_common >= 2:
        bootstrap_effects = bootstrap_numerator / total_weight
        permutation_effects = permutation_numerator / total_weight
        ci_lower, ci_upper = _percentile_interval(
            bootstrap_effects,
            confidence_level,
        )
        permutation_p_value = _two_sided_permutation_p_value(
            permutation_effects,
            effect,
        )
        valid_bootstrap_replicates = len(bootstrap_effects)
        valid_permutation_replicates = len(permutation_effects)

    return {
        "effect": effect,
        "standardized_group_mean": standardized_group_mean,
        "standardized_non_group_mean": standardized_non_group_mean,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper,
        "permutation_p_value": permutation_p_value,
        "group_n_cells": n_group_common,
        "non_group_n_cells": n_non_group_common,
        "group_retained_pct": group_retained_pct,
        "non_group_retained_pct": non_group_retained_pct,
        "mixed_bin_count": mixed_bin_count,
        "within_bin_coverage_difference_reads": within_bin_coverage_difference,
        "valid_bootstrap_replicates": valid_bootstrap_replicates,
        "valid_permutation_replicates": valid_permutation_replicates,
        "status": ";".join(status_parts) if status_parts else "ok",
    }


def _bootstrap_one_amplicon(job: Dict) -> Dict:
    """Compute estimates and paired bootstrap intervals for one amplicon."""
    amplicon = job["amplicon"]
    amplicon_index = job["amplicon_index"]
    mod_values = np.asarray(job["mod_values"], dtype=float)
    count_values = np.asarray(job["count_values"], dtype=float)
    hq_mask = np.asarray(job["hq_mask"], dtype=bool)
    min_reads = job["min_reads"]
    bootstrap_iterations = job["bootstrap_iterations"]
    permutation_iterations = job["permutation_iterations"]
    confidence_level = job["confidence_level"]
    base_seed = job["seed"]
    batch_size = job["batch_size"]
    coverage_exact_max_reads = job["coverage_exact_max_reads"]
    coverage_bin_width_reads = job["coverage_bin_width_reads"]

    valid_all = np.isfinite(mod_values) & np.isfinite(count_values) & (count_values >= min_reads)
    valid_hq = valid_all & hq_mask
    n_all = int(valid_all.sum())
    n_hq = int(valid_hq.sum())
    n_non_hq = n_all - n_hq
    all_estimate = _point_estimate(mod_values, valid_all)
    hq_estimate = _point_estimate(mod_values, valid_hq)
    non_hq_estimate = _point_estimate(mod_values, valid_all & ~hq_mask)
    delta_estimate = hq_estimate - all_estimate if np.isfinite(all_estimate) and np.isfinite(hq_estimate) else np.nan
    hq_non_hq_delta = (
        hq_estimate - non_hq_estimate
        if np.isfinite(hq_estimate) and np.isfinite(non_hq_estimate)
        else np.nan
    )

    all_values = mod_values[valid_all]
    hq_values = mod_values[valid_hq]
    non_hq_values = mod_values[valid_all & ~hq_mask]
    eligible_count_values = count_values[valid_all]
    eligible_hq_mask = hq_mask[valid_all]
    seed_sequence = np.random.SeedSequence([base_seed, amplicon_index])
    (
        all_seed,
        hq_seed,
        delta_seed,
        permutation_seed,
        coverage_bootstrap_seed,
        coverage_permutation_seed,
    ) = seed_sequence.spawn(6)

    all_bootstrap = np.array([], dtype=float)
    if n_all >= 2:
        all_bootstrap = _bootstrap_means(
            all_values, bootstrap_iterations, np.random.default_rng(all_seed), batch_size
        )

    hq_bootstrap = np.array([], dtype=float)
    if n_hq >= 2:
        hq_bootstrap = _bootstrap_means(
            hq_values, bootstrap_iterations, np.random.default_rng(hq_seed), batch_size
        )

    delta_bootstrap = np.array([], dtype=float)
    if n_all >= 2 and n_hq >= 2:
        delta_bootstrap = _bootstrap_stratified_delta(
            hq_values,
            non_hq_values,
            bootstrap_iterations,
            np.random.default_rng(delta_seed),
            batch_size,
        )

    permutation_null = np.array([], dtype=float)
    permutation_selected_means = np.array([], dtype=float)
    permutation_p_value = np.nan
    if n_hq >= 2 and n_non_hq >= 2:
        permutation_null = _permutation_hq_minus_all(
            all_values,
            n_hq,
            permutation_iterations,
            np.random.default_rng(permutation_seed),
            batch_size,
        )
        permutation_selected_means = permutation_null + all_estimate
        permutation_null = permutation_selected_means - all_estimate
        permutation_p_value = _two_sided_permutation_p_value(
            permutation_null, delta_estimate
        )

    all_lower, all_upper = _percentile_interval(all_bootstrap, confidence_level)
    hq_lower, hq_upper = _percentile_interval(hq_bootstrap, confidence_level)
    delta_lower, delta_upper = _percentile_interval(delta_bootstrap, confidence_level)
    hq_non_hq_bootstrap = np.array([], dtype=float)
    if len(delta_bootstrap) and n_non_hq >= 2:
        hq_non_hq_bootstrap = delta_bootstrap * n_all / float(n_non_hq)
    hq_non_hq_lower, hq_non_hq_upper = _percentile_interval(
        hq_non_hq_bootstrap,
        confidence_level,
    )

    coverage_adjusted = _coverage_adjusted_resampling(
        mod_values=all_values,
        count_values=eligible_count_values,
        group_mask=eligible_hq_mask,
        bootstrap_iterations=bootstrap_iterations,
        permutation_iterations=permutation_iterations,
        confidence_level=confidence_level,
        exact_max_reads=coverage_exact_max_reads,
        bin_width_reads=coverage_bin_width_reads,
        bootstrap_rng=np.random.default_rng(coverage_bootstrap_seed),
        permutation_rng=np.random.default_rng(coverage_permutation_seed),
    )

    status_parts = []
    if n_all < 2:
        status_parts.append("insufficient_all_cells")
    if n_hq < 2:
        status_parts.append("insufficient_in_group_cells")
    if n_non_hq < 2:
        status_parts.append("insufficient_out_group_cells")

    return {
        "amplicon": amplicon,
        "all_cells_estimate_pct": all_estimate,
        "all_cells_ci_lower_pct": all_lower,
        "all_cells_ci_upper_pct": all_upper,
        "all_cells_n_cells": n_all,
        "in_group_estimate_pct": hq_estimate,
        "in_group_ci_lower_pct": hq_lower,
        "in_group_ci_upper_pct": hq_upper,
        "in_group_n_cells": n_hq,
        "out_group_estimate_pct": non_hq_estimate,
        "out_group_n_cells": n_non_hq,
        "in_group_minus_all_cells_pct": delta_estimate,
        "in_group_minus_out_group_pct": hq_non_hq_delta,
        "in_group_minus_out_group_ci_lower_pct": hq_non_hq_lower,
        "in_group_minus_out_group_ci_upper_pct": hq_non_hq_upper,
        "in_group_minus_all_cells_ci_lower_pct": delta_lower,
        "in_group_minus_all_cells_ci_upper_pct": delta_upper,
        "permutation_p_value": permutation_p_value,
        "bh_adjusted_p_value": np.nan,
        "coverage_standardized_in_group_mean_pct": coverage_adjusted[
            "standardized_group_mean"
        ],
        "coverage_standardized_out_group_mean_pct": coverage_adjusted[
            "standardized_non_group_mean"
        ],
        "coverage_adjusted_in_group_minus_out_group_pct": coverage_adjusted["effect"],
        "coverage_adjusted_ci_lower_pct": coverage_adjusted["ci_lower"],
        "coverage_adjusted_ci_upper_pct": coverage_adjusted["ci_upper"],
        "coverage_adjusted_permutation_p_value": coverage_adjusted["permutation_p_value"],
        "coverage_adjusted_bh_p_value": np.nan,
        "coverage_adjusted_in_group_n_cells": coverage_adjusted["group_n_cells"],
        "coverage_adjusted_out_group_n_cells": coverage_adjusted["non_group_n_cells"],
        "coverage_adjusted_in_group_retained_pct": coverage_adjusted["group_retained_pct"],
        "coverage_adjusted_out_group_retained_pct": coverage_adjusted["non_group_retained_pct"],
        "coverage_adjusted_mixed_bin_count": coverage_adjusted["mixed_bin_count"],
        "coverage_adjusted_within_bin_coverage_difference_reads": coverage_adjusted[
            "within_bin_coverage_difference_reads"
        ],
        "valid_coverage_adjusted_bootstrap_replicates": coverage_adjusted[
            "valid_bootstrap_replicates"
        ],
        "valid_coverage_adjusted_permutation_replicates": coverage_adjusted[
            "valid_permutation_replicates"
        ],
        "coverage_exact_max_reads": coverage_exact_max_reads,
        "coverage_bin_width_reads": coverage_bin_width_reads,
        "coverage_adjusted_status": coverage_adjusted["status"],
        "valid_all_cells_bootstrap_replicates": len(all_bootstrap),
        "valid_in_group_bootstrap_replicates": len(hq_bootstrap),
        "valid_in_group_minus_all_cells_bootstrap_replicates": len(delta_bootstrap),
        "valid_permutation_replicates": len(permutation_null),
        "bootstrap_iterations": bootstrap_iterations,
        "permutation_iterations": permutation_iterations,
        "confidence_level": confidence_level,
        "seed": base_seed,
        "status": ";".join(status_parts) if status_parts else "ok",
        "_unconditional_permutation_selected_means": (
            permutation_selected_means
            if job.get("retain_permutation_distribution", False)
            else np.array([], dtype=float)
        ),
    }


def _compute_editing_rate_resampling_analyses(
    editing_summary: pd.DataFrame,
    quality_scores: pd.DataFrame,
    high_quality_codes: Iterable[str],
    min_reads_per_amplicon_per_cell: int,
    config: EditingRateCIConfig,
    n_processes: int = 1,
    include_permutation_outputs: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compute editing-rate intervals and their unconditional permutation nulls."""
    if "Color" not in quality_scores.columns:
        raise ValueError("Amplicon score table must contain a 'Color' column")
    if config.coverage_exact_max_reads < 0:
        raise ValueError("coverage_exact_max_reads must be >= 0")
    if config.coverage_bin_width_reads < 1:
        raise ValueError("coverage_bin_width_reads must be >= 1")

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
                "mod_values": categorical_mod_pct_to_cell_edit_pct(
                    pd.to_numeric(editing_summary[mod_column], errors="coerce").to_numpy(dtype=float),
                    amplicon=amplicon,
                    barcodes=editing_summary.index,
                ),
                "count_values": pd.to_numeric(editing_summary[count_column], errors="coerce").to_numpy(dtype=float),
                "hq_mask": hq_mask,
                "min_reads": min_reads_per_amplicon_per_cell,
                "bootstrap_iterations": config.bootstrap_iterations,
                "permutation_iterations": config.permutation_iterations,
                "confidence_level": config.confidence_level,
                "seed": config.seed,
                "batch_size": config.batch_size,
                "coverage_exact_max_reads": config.coverage_exact_max_reads,
                "coverage_bin_width_reads": config.coverage_bin_width_reads,
                "retain_permutation_distribution": include_permutation_outputs,
            }
        )

    if not jobs:
        return (
            pd.DataFrame(columns=RESULT_COLUMNS),
            pd.DataFrame(columns=UNCONDITIONAL_PERMUTATION_SUMMARY_COLUMNS),
            pd.DataFrame(columns=UNCONDITIONAL_PERMUTATION_SIMULATION_COLUMNS),
        )

    worker_count = min(max(1, int(n_processes)), len(jobs))
    if worker_count > 1:
        with mp.Pool(worker_count) as pool:
            results = pool.map(_bootstrap_one_amplicon, jobs)
    else:
        results = [_bootstrap_one_amplicon(job) for job in jobs]
    result_frame = pd.DataFrame(results, columns=RESULT_COLUMNS)
    result_frame["bh_adjusted_p_value"] = _benjamini_hochberg(
        result_frame["permutation_p_value"]
    )
    result_frame["coverage_adjusted_bh_p_value"] = _benjamini_hochberg(
        result_frame["coverage_adjusted_permutation_p_value"]
    )
    if not include_permutation_outputs:
        return (
            result_frame,
            pd.DataFrame(columns=UNCONDITIONAL_PERMUTATION_SUMMARY_COLUMNS),
            pd.DataFrame(columns=UNCONDITIONAL_PERMUTATION_SIMULATION_COLUMNS),
        )

    summary_rows = []
    simulation_rows = []
    for result, (_, row) in zip(results, result_frame.iterrows()):
        permutation_means = np.asarray(
            result["_unconditional_permutation_selected_means"], dtype=float
        )
        n_all = int(row["all_cells_n_cells"])
        n_selected = int(row["in_group_n_cells"])
        n_non_selected = int(row["out_group_n_cells"])
        all_estimate = float(row["all_cells_estimate_pct"])
        selected_estimate = float(row["in_group_estimate_pct"])
        status_parts = []
        if n_all < 2:
            status_parts.append("insufficient_all_cells")
        if n_selected < 2:
            status_parts.append("insufficient_selected_group_cells")
        if n_non_selected < 2:
            status_parts.append("insufficient_non_selected_cells")
        if len(permutation_means) and np.all(permutation_means == permutation_means[0]):
            status_parts.append("invariant_permutation_distribution")

        lower, upper = _percentile_interval(
            permutation_means, config.confidence_level
        )
        extreme_count = 0
        percentile = np.nan
        if len(permutation_means):
            observed_distance = abs(selected_estimate - all_estimate)
            extreme_count = int(
                np.count_nonzero(
                    np.abs(permutation_means - all_estimate) >= observed_distance
                )
            )
            percentile = float(
                100.0 * np.mean(permutation_means <= selected_estimate)
            )

        summary_rows.append(
            {
                "amplicon": row["amplicon"],
                "all_estimate_pct": all_estimate,
                "all_n_cells": n_all,
                "selected_group_estimate_pct": selected_estimate,
                "selected_group_n_cells": n_selected,
                "selected_group_minus_all_pct": row["in_group_minus_all_cells_pct"],
                "permuted_median_pct": (
                    float(np.median(permutation_means))
                    if len(permutation_means)
                    else np.nan
                ),
                "permuted_ci_lower_pct": lower,
                "permuted_ci_upper_pct": upper,
                "selected_group_permutation_percentile": percentile,
                "permutation_p_value": row["permutation_p_value"],
                "bh_adjusted_p_value": row["bh_adjusted_p_value"],
                "extreme_permutation_count": extreme_count,
                "valid_permutations": len(permutation_means),
                "requested_permutations": config.permutation_iterations,
                "confidence_level": config.confidence_level,
                "seed": config.seed,
                "status": ";".join(status_parts) if status_parts else "ok",
            }
        )
        simulation_rows.extend(
            {
                "amplicon": row["amplicon"],
                "permutation_index": permutation_index + 1,
                "permuted_selected_estimate_pct": permuted_estimate,
                "selected_group_n_cells": n_selected,
                "eligible_all_n_cells": n_all,
                "seed": config.seed,
            }
            for permutation_index, permuted_estimate in enumerate(permutation_means)
        )

    summaries = pd.DataFrame(summary_rows).reindex(
        columns=UNCONDITIONAL_PERMUTATION_SUMMARY_COLUMNS
    )
    simulations = pd.DataFrame(simulation_rows).reindex(
        columns=UNCONDITIONAL_PERMUTATION_SIMULATION_COLUMNS
    )
    return result_frame, summaries, simulations


def compute_editing_rate_confidence_intervals(
    editing_summary: pd.DataFrame,
    quality_scores: pd.DataFrame,
    high_quality_codes: Iterable[str],
    min_reads_per_amplicon_per_cell: int,
    config: EditingRateCIConfig,
    n_processes: int = 1,
) -> pd.DataFrame:
    """Compute pointwise editing-rate intervals for all and the configured group."""
    results, _, _ = _compute_editing_rate_resampling_analyses(
        editing_summary,
        quality_scores,
        high_quality_codes,
        min_reads_per_amplicon_per_cell,
        config,
        n_processes=n_processes,
        include_permutation_outputs=False,
    )
    return results


def compute_editing_rate_resampling_analyses(
    editing_summary: pd.DataFrame,
    quality_scores: pd.DataFrame,
    high_quality_codes: Iterable[str],
    min_reads_per_amplicon_per_cell: int,
    config: EditingRateCIConfig,
    n_processes: int = 1,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compute intervals plus the stored unconditional permutation distribution."""
    return _compute_editing_rate_resampling_analyses(
        editing_summary,
        quality_scores,
        high_quality_codes,
        min_reads_per_amplicon_per_cell,
        config,
        n_processes=n_processes,
        include_permutation_outputs=True,
    )
