"""Bootstrap confidence intervals and significance tests for editing rates."""

from dataclasses import dataclass
import logging
import math
import multiprocessing as mp
import os
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


RESULT_COLUMNS = [
    "amplicon",
    "all_cells_estimate_pct",
    "all_cells_ci_lower_pct",
    "all_cells_ci_upper_pct",
    "all_cells_n_cells",
    "in_group_estimate_pct",
    "in_group_ci_lower_pct",
    "in_group_ci_upper_pct",
    "in_group_n_cells",
    "out_group_estimate_pct",
    "out_group_n_cells",
    "in_group_minus_all_cells_pct",
    "in_group_minus_out_group_pct",
    "in_group_minus_out_group_ci_lower_pct",
    "in_group_minus_out_group_ci_upper_pct",
    "in_group_minus_all_cells_ci_lower_pct",
    "in_group_minus_all_cells_ci_upper_pct",
    "permutation_p_value",
    "bh_adjusted_p_value",
    "coverage_standardized_in_group_mean_pct",
    "coverage_standardized_out_group_mean_pct",
    "coverage_adjusted_in_group_minus_out_group_pct",
    "coverage_adjusted_ci_lower_pct",
    "coverage_adjusted_ci_upper_pct",
    "coverage_adjusted_permutation_p_value",
    "coverage_adjusted_bh_p_value",
    "coverage_adjusted_in_group_n_cells",
    "coverage_adjusted_out_group_n_cells",
    "coverage_adjusted_in_group_retained_pct",
    "coverage_adjusted_out_group_retained_pct",
    "coverage_adjusted_mixed_bin_count",
    "coverage_adjusted_within_bin_coverage_difference_reads",
    "valid_coverage_adjusted_bootstrap_replicates",
    "valid_coverage_adjusted_permutation_replicates",
    "coverage_exact_max_reads",
    "coverage_bin_width_reads",
    "coverage_adjusted_status",
    "valid_all_cells_bootstrap_replicates",
    "valid_in_group_bootstrap_replicates",
    "valid_in_group_minus_all_cells_bootstrap_replicates",
    "valid_permutation_replicates",
    "bootstrap_iterations",
    "permutation_iterations",
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

FIXED_CELL_COUNT_DEPTHS = (100, 200, 400, 800, 1000)

FIXED_CELL_DEPTH_STABILITY_COLUMNS = [
    "amplicon",
    "cohort",
    "requested_sample_n_cells",
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
    "status",
]

UNCONDITIONAL_PERMUTATION_SUMMARY_COLUMNS = [
    "amplicon",
    "all_estimate_pct",
    "all_n_cells",
    "selected_group_estimate_pct",
    "selected_group_n_cells",
    "selected_group_minus_all_pct",
    "permuted_median_pct",
    "permuted_ci_lower_pct",
    "permuted_ci_upper_pct",
    "selected_group_permutation_percentile",
    "permutation_p_value",
    "bh_adjusted_p_value",
    "extreme_permutation_count",
    "valid_permutations",
    "requested_permutations",
    "confidence_level",
    "seed",
    "status",
]

UNCONDITIONAL_PERMUTATION_SIMULATION_COLUMNS = [
    "amplicon",
    "permutation_index",
    "permuted_selected_estimate_pct",
    "selected_group_n_cells",
    "eligible_all_n_cells",
    "seed",
]


CI_PLOT_SUFFIXES = (
    ".10_EditingRateConfidenceIntervals",
    ".11_EditingRateQualityDelta",
    ".14_EditingRateCoverageAdjustedEffects",
    ".16_EditingRateUnconditionalPermutation",
)

DEPTH_STABILITY_PLOT_SUFFIXES = (
    ".12_EditingRateDepthStability",
    ".13_EditingRateRelativeDepthStability",
    ".15_EditingRateFixedCellDepthStability",
)


@dataclass(frozen=True)
class EditingRateCIConfig:
    """Configuration for editing-rate bootstrap confidence intervals."""

    enabled: bool = True
    bootstrap_iterations: int = 10_000
    permutation_iterations: int = 10_000
    confidence_level: float = 0.95
    seed: int = 42
    batch_size: int = 64
    coverage_exact_max_reads: int = 10
    coverage_bin_width_reads: int = 5


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
                "mod_values": pd.to_numeric(editing_summary[mod_column], errors="coerce").to_numpy(dtype=float),
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


def _empty_fixed_cell_depth_stability_row(
    amplicon: str,
    cohort: str,
    requested_sample_n_cells: int,
    eligible_n_cells: int,
    full_estimate: float,
    job: Dict,
) -> Dict:
    """Describe a fixed cell-count request that cannot be sampled exactly."""
    return {
        "amplicon": amplicon,
        "cohort": cohort,
        "requested_sample_n_cells": requested_sample_n_cells,
        "sample_n_cells": 0,
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
        "status": "insufficient_cells",
    }


def _fixed_cell_depth_stability_one_amplicon(job: Dict) -> List[Dict]:
    """Evaluate exact fixed cell-count subsamples for one amplicon."""
    amplicon = job["amplicon"]
    mod_values = np.asarray(job["mod_values"], dtype=float)
    count_values = np.asarray(job["count_values"], dtype=float)
    hq_mask = np.asarray(job["hq_mask"], dtype=bool)
    valid_all = (
        np.isfinite(mod_values)
        & np.isfinite(count_values)
        & (count_values >= job["min_reads"])
    )
    requested_sizes = tuple(int(value) for value in job["fixed_sample_sizes"])
    amplicon_seed = np.random.SeedSequence([job["seed"], job["amplicon_index"], 1])
    cohort_seeds = amplicon_seed.spawn(2)
    results: List[Dict] = []

    for cohort_index, (cohort, cohort_mask) in enumerate(
        (("all", valid_all), ("hq", valid_all & hq_mask))
    ):
        values = mod_values[cohort_mask]
        eligible_n_cells = len(values)
        full_estimate = float(np.mean(values)) if eligible_n_cells else np.nan
        available_sizes = [
            sample_size
            for sample_size in requested_sizes
            if eligible_n_cells >= sample_size and eligible_n_cells >= 2
        ]
        sampled_means = (
            _nested_subsample_means(
                values,
                available_sizes,
                job["iterations"],
                np.random.default_rng(cohort_seeds[cohort_index]),
            )
            if available_sizes
            else {}
        )

        for requested_sample_n_cells in requested_sizes:
            if requested_sample_n_cells not in sampled_means:
                results.append(
                    _empty_fixed_cell_depth_stability_row(
                        amplicon,
                        cohort,
                        requested_sample_n_cells,
                        eligible_n_cells,
                        full_estimate,
                        job,
                    )
                )
                continue

            means = sampled_means[requested_sample_n_cells]
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
                    "requested_sample_n_cells": requested_sample_n_cells,
                    "sample_n_cells": requested_sample_n_cells,
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
                    "status": "ok",
                }
            )
    return results


def compute_editing_rate_fixed_cell_depth_stability(
    editing_summary: pd.DataFrame,
    quality_scores: pd.DataFrame,
    high_quality_codes: Iterable[str],
    min_reads_per_amplicon_per_cell: int,
    config: EditingRateDepthStabilityConfig,
    n_processes: int = 1,
) -> pd.DataFrame:
    """Measure editing-rate stability at exact fixed eligible-cell counts."""
    if "Color" not in quality_scores.columns:
        raise ValueError("Amplicon score table must contain a 'Color' column")

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
                "fixed_sample_sizes": FIXED_CELL_COUNT_DEPTHS,
                "confidence_level": config.confidence_level,
                "seed": config.seed,
            }
        )

    if not jobs:
        return pd.DataFrame(columns=FIXED_CELL_DEPTH_STABILITY_COLUMNS)

    worker_count = min(max(1, int(n_processes)), len(jobs))
    if worker_count > 1:
        with mp.Pool(worker_count) as pool:
            nested_results = pool.map(_fixed_cell_depth_stability_one_amplicon, jobs)
    else:
        nested_results = [_fixed_cell_depth_stability_one_amplicon(job) for job in jobs]
    results = [row for amplicon_rows in nested_results for row in amplicon_rows]
    return pd.DataFrame(results).reindex(columns=FIXED_CELL_DEPTH_STABILITY_COLUMNS)


def _finite_interval_rows(results: pd.DataFrame, prefix: str) -> pd.Series:
    return (
        pd.to_numeric(results[f"{prefix}_estimate_pct"], errors="coerce").notna()
        & pd.to_numeric(results[f"{prefix}_ci_lower_pct"], errors="coerce").notna()
        & pd.to_numeric(results[f"{prefix}_ci_upper_pct"], errors="coerce").notna()
    )


def _significant_plot_rows(results: pd.DataFrame) -> pd.DataFrame:
    """Return amplicons passing the coverage-adjusted BH threshold."""
    adjusted_p_values = pd.to_numeric(
        results["coverage_adjusted_bh_p_value"],
        errors="coerce",
    )
    significant = adjusted_p_values.notna() & (adjusted_p_values <= 0.05)
    return results.loc[significant].copy()


def _remove_plot_artifacts(output_root: str, suffixes: Sequence[str]) -> None:
    """Remove PNG/PDF artifacts that may have been created by an earlier run."""
    for suffix in suffixes:
        for extension in (".png", ".pdf"):
            path = output_root + suffix + extension
            try:
                os.remove(path)
            except FileNotFoundError:
                pass


def remove_editing_rate_plot_artifacts(output_root: str) -> None:
    """Remove every generated editing-rate plot from an earlier run."""
    _remove_plot_artifacts(
        output_root,
        CI_PLOT_SUFFIXES + DEPTH_STABILITY_PLOT_SUFFIXES,
    )


def _row_plot_height(n_rows: int, *, legend: bool = False) -> float:
    """Scale a horizontal interval plot to its number of amplicons."""
    overhead = 2.4 if legend else 2.0
    minimum = 3.8 if legend else 3.2
    return max(minimum, overhead + (0.4 * max(1, n_rows)))


def _interval_axis_limits(
    lower_values: Sequence[float],
    upper_values: Sequence[float],
    *,
    minimum_span: float,
    include_zero: bool = False,
    bounds: Optional[Tuple[float, float]] = None,
) -> Optional[Tuple[float, float]]:
    """Return padded, readable limits around finite confidence intervals."""
    lower = np.asarray(lower_values, dtype=float)
    upper = np.asarray(upper_values, dtype=float)
    finite_values = np.concatenate(
        [lower[np.isfinite(lower)], upper[np.isfinite(upper)]]
    )
    if finite_values.size == 0:
        return None

    data_lower = float(np.min(finite_values))
    data_upper = float(np.max(finite_values))
    if include_zero:
        data_lower = min(0.0, data_lower)
        data_upper = max(0.0, data_upper)

    data_span = data_upper - data_lower
    padded_span = max(minimum_span, data_span * 1.16)
    center = (data_lower + data_upper) / 2.0
    axis_lower = center - (padded_span / 2.0)
    axis_upper = center + (padded_span / 2.0)

    if bounds is not None:
        bound_lower, bound_upper = bounds
        bounded_span = bound_upper - bound_lower
        if padded_span >= bounded_span:
            return float(bound_lower), float(bound_upper)
        if axis_lower < bound_lower:
            axis_upper += bound_lower - axis_lower
            axis_lower = bound_lower
        if axis_upper > bound_upper:
            axis_lower -= axis_upper - bound_upper
            axis_upper = bound_upper
        axis_lower = max(bound_lower, axis_lower)
        axis_upper = min(bound_upper, axis_upper)

    return float(axis_lower), float(axis_upper)


def _depth_stability_layout(n_amplicons: int) -> Tuple[int, int, float, float]:
    """Choose a compact facet grid and figure size for a stability plot."""
    if n_amplicons < 1:
        raise ValueError("At least one amplicon is required for a stability plot")
    if n_amplicons <= 3:
        n_columns = n_amplicons
    elif n_amplicons == 4:
        n_columns = 2
    else:
        n_columns = 3
    n_rows = int(math.ceil(n_amplicons / n_columns))
    figure_width = max(8.0, 6.0 * n_columns)
    figure_height = 1.5 + (3.8 * n_rows)
    return n_rows, n_columns, figure_width, figure_height


def _write_coverage_adjusted_effect_plot(
    results: pd.DataFrame,
    output_root: str,
) -> List[Dict[str, str]]:
    """Compare raw and coverage-adjusted InGroup-minus-OutGroup effects."""
    required_columns = [
        "in_group_minus_out_group_pct",
        "in_group_minus_out_group_ci_lower_pct",
        "in_group_minus_out_group_ci_upper_pct",
        "coverage_adjusted_in_group_minus_out_group_pct",
        "coverage_adjusted_ci_lower_pct",
        "coverage_adjusted_ci_upper_pct",
        "bh_adjusted_p_value",
        "coverage_adjusted_bh_p_value",
    ]
    plotted = results.copy()
    for column in required_columns:
        plotted[column] = pd.to_numeric(plotted[column], errors="coerce")
    raw_valid = (
        plotted["in_group_minus_out_group_pct"].notna()
        & plotted["in_group_minus_out_group_ci_lower_pct"].notna()
        & plotted["in_group_minus_out_group_ci_upper_pct"].notna()
    )
    adjusted_valid = (
        plotted["coverage_adjusted_in_group_minus_out_group_pct"].notna()
        & plotted["coverage_adjusted_ci_lower_pct"].notna()
        & plotted["coverage_adjusted_ci_upper_pct"].notna()
    )
    plotted = plotted.loc[raw_valid | adjusted_valid].copy()
    if plotted.empty:
        return []

    plotted = plotted.sort_values(
        "coverage_adjusted_in_group_minus_out_group_pct",
        ascending=True,
        na_position="first",
    ).reset_index(drop=True)
    fig_height = _row_plot_height(len(plotted), legend=True)
    fig, ax = plt.subplots(figsize=(12, fig_height))
    y_positions = np.arange(len(plotted), dtype=float)
    series = [
        (
            "in_group_minus_out_group_pct",
            "in_group_minus_out_group_ci_lower_pct",
            "in_group_minus_out_group_ci_upper_pct",
            "bh_adjusted_p_value",
            "Gray circle: InGroup minus OutGroup (unconditional BH)",
            "#7F7F7F",
            "o",
            -0.12,
        ),
        (
            "coverage_adjusted_in_group_minus_out_group_pct",
            "coverage_adjusted_ci_lower_pct",
            "coverage_adjusted_ci_upper_pct",
            "coverage_adjusted_bh_p_value",
            "Blue square: coverage-adjusted InGroup minus OutGroup (coverage-adjusted BH)",
            "#4C78A8",
            "s",
            0.12,
        ),
    ]
    legend_handles = []
    for estimate_column, lower_column, upper_column, p_column, label, color, marker, offset in series:
        valid = (
            plotted[estimate_column].notna()
            & plotted[lower_column].notna()
            & plotted[upper_column].notna()
        )
        for row_index in np.flatnonzero(valid.to_numpy()):
            row = plotted.iloc[row_index]
            estimate = float(row[estimate_column])
            lower = float(row[lower_column])
            upper = float(row[upper_column])
            significant = np.isfinite(row[p_column]) and row[p_column] <= 0.05
            ax.errorbar(
                estimate,
                y_positions[row_index] + offset,
                xerr=np.array([[max(0.0, estimate - lower)], [max(0.0, upper - estimate)]]),
                fmt=marker,
                capsize=3,
                color=color,
                markerfacecolor=color if significant else "white",
                markeredgecolor=color,
            )
        if valid.any():
            # Use a separate, always-open handle so the legend does not inherit
            # the significance state of whichever amplicon happened to plot first.
            legend_handle = ax.errorbar(
                [],
                [],
                xerr=np.empty((2, 0)),
                fmt=marker,
                capsize=3,
                color=color,
                markerfacecolor="white",
                markeredgecolor=color,
            )
            legend_handles.append((legend_handle, label))

    ax.axvline(0, color="black", linestyle="--", linewidth=1)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(plotted["amplicon"])
    ax.set_xlabel("InGroup minus OutGroup (percentage points)")
    ax.set_ylabel("Amplicon")
    ax.set_title("Raw and Coverage-Adjusted Editing-Rate Effects")
    if legend_handles:
        ax.legend(
            [item[0] for item in legend_handles],
            [item[1] for item in legend_handles],
            title="Open = not significant; filled = corresponding BH-adjusted p-value ≤ 0.05",
        )
    interval_limits = _interval_axis_limits(
        np.concatenate(
            [
                plotted["in_group_minus_out_group_ci_lower_pct"].to_numpy(dtype=float),
                plotted["coverage_adjusted_ci_lower_pct"].to_numpy(dtype=float),
            ]
        ),
        np.concatenate(
            [
                plotted["in_group_minus_out_group_ci_upper_pct"].to_numpy(dtype=float),
                plotted["coverage_adjusted_ci_upper_pct"].to_numpy(dtype=float),
            ]
        ),
        minimum_span=2.0,
        include_zero=True,
    )
    if interval_limits is not None:
        ax.set_xlim(*interval_limits)
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    plot_root = output_root + ".14_EditingRateCoverageAdjustedEffects"
    fig.savefig(plot_root + ".pdf", bbox_inches="tight")
    fig.savefig(plot_root + ".png", bbox_inches="tight")
    plt.close(fig)
    return [
        {
            "plot_name": plot_root,
            "plot_title": "Raw and coverage-adjusted editing-rate effects",
            "plot_label": (
                "Bootstrap intervals for the unconditional and coverage-stratified "
                "configured-group-minus-remaining-cell effects. Filled markers "
                "denote significance after separate BH adjustments across all "
                "testable amplicons."
            ),
        }
    ]


def write_editing_rate_unconditional_permutation_plot(
    summaries: pd.DataFrame,
    simulations: pd.DataFrame,
    output_root: str,
) -> List[Dict[str, str]]:
    """Plot unconditional permutation distributions with selected-group estimates."""
    suffix = ".16_EditingRateUnconditionalPermutation"
    _remove_plot_artifacts(output_root, (suffix,))
    if summaries.empty or simulations.empty:
        return []

    plotted = summaries.copy()
    for column in (
        "all_estimate_pct",
        "selected_group_estimate_pct",
        "valid_permutations",
    ):
        plotted[column] = pd.to_numeric(plotted[column], errors="coerce")
    plotted = plotted.loc[
        plotted["all_estimate_pct"].notna()
        & plotted["selected_group_estimate_pct"].notna()
        & (plotted["valid_permutations"] > 0)
    ].copy()
    if plotted.empty:
        return []
    plotted = plotted.sort_values(
        ["all_estimate_pct", "amplicon"], ascending=[True, True]
    ).reset_index(drop=True)
    permutation_groups = {
        str(amplicon): pd.to_numeric(
            group["permuted_selected_estimate_pct"], errors="coerce"
        )
        .dropna()
        .to_numpy(dtype=float)
        for amplicon, group in simulations.groupby("amplicon", sort=False)
    }
    plotted = plotted.loc[
        plotted["amplicon"].astype(str).map(
            lambda amplicon: len(permutation_groups.get(amplicon, ())) > 0
        )
    ].reset_index(drop=True)
    if plotted.empty:
        return []

    boxes = [permutation_groups[str(amplicon)] for amplicon in plotted["amplicon"]]
    selected_estimates = plotted["selected_group_estimate_pct"].to_numpy(dtype=float)
    all_plot_values = np.concatenate([np.concatenate(boxes), selected_estimates])
    figure_width = max(10.0, 3.0 + (0.55 * len(plotted)))
    fig, ax = plt.subplots(figsize=(figure_width, 8.0))
    boxplot = ax.boxplot(
        boxes,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#1F1F1F", "linewidth": 1.4},
        whiskerprops={"color": "#4C78A8"},
        capprops={"color": "#4C78A8"},
    )
    for box in boxplot["boxes"]:
        box.set(facecolor="#A6C8E0", edgecolor="#4C78A8", alpha=0.9)
    positions = np.arange(1, len(plotted) + 1, dtype=float)
    ax.scatter(
        positions,
        selected_estimates,
        marker="D",
        s=44,
        color="#F58518",
        edgecolor="#1F1F1F",
        linewidth=0.6,
        zorder=4,
        label="Configured analysis-group estimate",
    )
    ax.set_xticks(positions)
    ax.set_xticklabels(plotted["amplicon"], rotation=65, ha="right", fontsize=9)
    limits = _interval_axis_limits(
        all_plot_values,
        all_plot_values,
        minimum_span=5.0,
        bounds=(0.0, 100.0),
    )
    if limits is not None:
        ax.set_ylim(*limits)
    ax.set_xlabel("Amplicon")
    ax.set_ylabel("Mean inferred allele editing percentage")
    ax.set_title("Unconditional Selected-Group Permutation Distribution")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    plot_root = output_root + suffix
    fig.savefig(plot_root + ".pdf", bbox_inches="tight")
    fig.savefig(plot_root + ".png", bbox_inches="tight")
    plt.close(fig)
    return [
        {
            "plot_name": plot_root,
            "plot_title": "Unconditional editing-rate permutation distribution",
            "plot_label": (
                "For each amplicon, boxplots show configured-group-sized subsets "
                "drawn without replacement from all eligible cells. Orange diamonds "
                "show the observed configured analysis-group estimates."
            ),
        }
    ]


def write_editing_rate_ci_plots(results: pd.DataFrame, output_root: str) -> List[Dict[str, str]]:
    """Plot controlled-significant intervals and all adjusted-effect comparisons."""
    plot_metadata: List[Dict[str, str]] = []
    _remove_plot_artifacts(output_root, CI_PLOT_SUFFIXES)
    if results.empty:
        return plot_metadata

    comparison_metadata = _write_coverage_adjusted_effect_plot(results, output_root)

    ordered = _significant_plot_rows(results)
    if ordered.empty:
        logging.info(
            "No amplicons passed the coverage-adjusted BH threshold of 0.05; "
            "skipping significance-filtered editing-rate confidence interval plots"
        )
        return comparison_metadata

    ordered["all_cells_estimate_pct"] = pd.to_numeric(
        ordered["all_cells_estimate_pct"], errors="coerce"
    )
    ordered = ordered.sort_values(
        "all_cells_estimate_pct", ascending=True, na_position="first"
    ).reset_index(drop=True)

    all_valid = _finite_interval_rows(ordered, "all_cells")
    hq_valid = _finite_interval_rows(ordered, "in_group")
    if all_valid.any() or hq_valid.any():
        fig_height = _row_plot_height(len(ordered), legend=True)
        fig, ax = plt.subplots(figsize=(12, fig_height))
        y_positions = np.arange(len(ordered), dtype=float)
        for valid, prefix, label, color, offset in [
            (all_valid, "all_cells", "AllCells", "#4C78A8", -0.12),
            (hq_valid, "in_group", "InGroup", "#F58518", 0.12),
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
        interval_limits = _interval_axis_limits(
            np.concatenate(
                [
                    ordered.loc[all_valid, "all_cells_ci_lower_pct"].to_numpy(dtype=float),
                    ordered.loc[hq_valid, "in_group_ci_lower_pct"].to_numpy(dtype=float),
                ]
            ),
            np.concatenate(
                [
                    ordered.loc[all_valid, "all_cells_ci_upper_pct"].to_numpy(dtype=float),
                    ordered.loc[hq_valid, "in_group_ci_upper_pct"].to_numpy(dtype=float),
                ]
            ),
            minimum_span=5.0,
            bounds=(0.0, 100.0),
        )
        if interval_limits is not None:
            ax.set_xlim(*interval_limits)
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
                    "and configured analysis-group cells among amplicons with a "
                    "coverage-adjusted BH p-value at or below 0.05."
                ),
            }
        )

    delta_valid = (
        pd.to_numeric(ordered["in_group_minus_all_cells_pct"], errors="coerce").notna()
        & pd.to_numeric(
            ordered["in_group_minus_all_cells_ci_lower_pct"], errors="coerce"
        ).notna()
        & pd.to_numeric(
            ordered["in_group_minus_all_cells_ci_upper_pct"], errors="coerce"
        ).notna()
    )
    if delta_valid.any():
        fig_height = _row_plot_height(len(ordered))
        fig, ax = plt.subplots(figsize=(12, fig_height))
        subset = ordered.loc[delta_valid]
        positions = np.arange(len(ordered), dtype=float)[delta_valid.to_numpy()]
        estimate = subset["in_group_minus_all_cells_pct"].to_numpy(dtype=float)
        lower = subset["in_group_minus_all_cells_ci_lower_pct"].to_numpy(dtype=float)
        upper = subset["in_group_minus_all_cells_ci_upper_pct"].to_numpy(dtype=float)
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
        interval_limits = _interval_axis_limits(
            lower,
            upper,
            minimum_span=2.0,
            include_zero=True,
        )
        if interval_limits is not None:
            ax.set_xlim(*interval_limits)
        ax.set_xlabel("InGroup minus AllCells editing percentage points")
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
                    "Paired bootstrap interval for the configured analysis-group "
                    "estimate minus the all-analyzable-cell estimate among "
                    "amplicons with a coverage-adjusted BH p-value at or below "
                    "0.05."
                ),
            }
        )
    return plot_metadata + comparison_metadata


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
    x_column: str = "sample_percent",
    x_label: str = "Eligible cohort retained (%)",
    x_limits: Tuple[float, float] = (0.0, 102.0),
    x_tick_values: Optional[Sequence[float]] = None,
    annotate_missing: bool = False,
):
    """Create a stability facet grid scaled to the included amplicons."""
    n_rows, n_columns, figure_width, figure_height = _depth_stability_layout(
        len(amplicon_order)
    )
    if subtitle:
        figure_height += 0.7
    fig, axes = plt.subplots(
        n_rows,
        n_columns,
        figsize=(figure_width, figure_height),
        sharex=True,
        sharey=True,
        squeeze=False,
    )
    axes_flat = axes.ravel()
    finite_bounds = plotted[[lower_column, upper_column]].to_numpy(dtype=float)
    finite_bounds = finite_bounds[np.isfinite(finite_bounds)]
    maximum_absolute_bound = (
        float(np.max(np.abs(finite_bounds))) if finite_bounds.size else 0.0
    )
    y_limit = max(0.1, maximum_absolute_bound * 1.1)
    styles = {
        "all": ("All analyzable cells", "#4C78A8"),
        "hq": ("Configured analysis-group cells", "#F58518"),
    }
    legend_handles = []
    legend_labels = []

    for facet_index, (axis, amplicon) in enumerate(zip(axes_flat, amplicon_order)):
        amplicon_rows = plotted.loc[plotted["amplicon"] == amplicon]
        plotted_cohort = False
        for cohort, (label, color) in styles.items():
            cohort_rows = amplicon_rows.loc[
                amplicon_rows["cohort"] == cohort
            ].sort_values(x_column)
            if cohort_rows.empty:
                continue
            x_values = cohort_rows[x_column].to_numpy(dtype=float)
            medians = cohort_rows[median_column].to_numpy(dtype=float)
            lower = cohort_rows[lower_column].to_numpy(dtype=float)
            upper = cohort_rows[upper_column].to_numpy(dtype=float)
            line = axis.plot(
                x_values,
                medians,
                marker="o",
                markersize=3,
                linewidth=1.3,
                color=color,
                label=label,
            )[0]
            axis.fill_between(x_values, lower, upper, color=color, alpha=0.2)
            plotted_cohort = True
            if label not in legend_labels:
                legend_handles.append(line)
                legend_labels.append(label)
        if annotate_missing and not plotted_cohort:
            axis.text(
                0.5,
                0.5,
                "No cohort has at least 100 eligible cells",
                ha="center",
                va="center",
                transform=axis.transAxes,
                fontsize=8,
                color="#555555",
            )
        axis.axhline(0, color="black", linestyle="--", linewidth=0.8, alpha=0.7)
        axis.set_title(amplicon, fontsize=8)
        axis.set_xlim(*x_limits)
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

    tick_values = (
        list(x_tick_values)
        if x_tick_values is not None
        else sorted(plotted[x_column].dropna().unique())
    )
    for axis in axes_flat[:len(amplicon_order)]:
        axis.set_xticks(tick_values)

    for row_index in range(n_rows):
        row_start = row_index * n_columns
        populated_count = min(n_columns, len(amplicon_order) - row_start)
        if populated_count <= 0:
            continue
        label_column = 1 if populated_count >= 2 else 0
        axes[row_index, label_column].set_xlabel(x_label, fontsize=9, labelpad=5)
        axes[row_index, 0].set_ylabel(y_label, fontsize=8, labelpad=7)

    fig.suptitle(title, fontsize=16, y=0.99)
    if subtitle:
        fig.text(0.5, 0.935, subtitle, ha="center", va="top", fontsize=10)
        legend_y = 0.875
        layout_top = 0.78
    else:
        legend_y = 0.935
        layout_top = 0.85
    if legend_handles:
        fig.legend(
            legend_handles,
            legend_labels,
            loc="upper center",
            bbox_to_anchor=(0.5, legend_y),
            ncol=2,
        )
    fig.tight_layout(rect=(0.02, 0.01, 1.0, layout_top), h_pad=2.0, w_pad=1.0)
    return fig


def _save_depth_stability_figure(fig, plot_root: str) -> None:
    fig.savefig(plot_root + ".pdf", bbox_inches="tight")
    fig.savefig(plot_root + ".png", bbox_inches="tight")
    plt.close(fig)


def write_editing_rate_depth_stability_plot(
    results: pd.DataFrame,
    output_root: str,
    significant_amplicons: Optional[Iterable[str]] = None,
) -> List[Dict[str, str]]:
    """Write absolute and relative stability plots with optional filtering."""
    _remove_plot_artifacts(output_root, DEPTH_STABILITY_PLOT_SUFFIXES)
    if results.empty:
        return []

    significance_filter_applied = significant_amplicons is not None
    plotted = results.loc[
        results["status"].isin(["ok", "full_reference"])
    ].copy()
    if significance_filter_applied:
        significant_amplicons = {str(amplicon) for amplicon in significant_amplicons}
        plotted = plotted.loc[
            plotted["amplicon"].astype(str).isin(significant_amplicons)
        ].copy()
        if plotted.empty:
            logging.info(
                "No significant amplicons had usable depth-stability results; "
                "skipping editing-rate depth-stability plots"
            )
            return []
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
    significance_clause = (
        " among amplicons with a coverage-adjusted BH p-value at or below 0.05"
        if significance_filter_applied
        else ""
    )
    plot_metadata = [
        {
            "plot_name": absolute_plot_root,
            "plot_title": "Editing-rate cell-depth stability",
            "plot_label": (
                "Finite-cohort downsampling stability bands for all analyzable "
                f"and configured analysis-group cells{significance_clause}. Bands "
                "show "
                "sensitivity to retained cell depth within this run, not "
                "uncertainty across biological replicates."
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
            "Amplicons included when the full configured-group editing rate is "
            f"at least {threshold:g}%"
        ),
        y_label="Relative deviation from full-cohort editing rate (%)",
    )
    relative_plot_root = output_root + ".13_EditingRateRelativeDepthStability"
    _save_depth_stability_figure(relative_figure, relative_plot_root)
    relative_significance = "significant " if significance_filter_applied else ""
    plot_metadata.append(
        {
            "plot_name": relative_plot_root,
            "plot_title": "Relative editing-rate cell-depth stability",
            "plot_label": (
                "Relative finite-cohort downsampling stability bands for "
                f"{relative_significance}amplicons meeting the configured analysis-group "
                "editing-rate "
                "threshold. Each cohort is normalized to its own full-cohort "
                "editing rate."
            ),
        }
    )
    return plot_metadata


def _fixed_cell_depth_stability_amplicon_order(results: pd.DataFrame) -> List[str]:
    """Order fixed-count panels by configured-group then all-cell editing rate."""
    references = results[
        ["amplicon", "cohort", "full_estimate_pct"]
    ].drop_duplicates(["amplicon", "cohort"])
    references = references.copy()
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
        hq_estimate = (
            reference_lookup.at[amplicon, "hq"]
            if amplicon in reference_lookup.index and "hq" in reference_lookup
            else np.nan
        )
        all_estimate = (
            reference_lookup.at[amplicon, "all"]
            if amplicon in reference_lookup.index and "all" in reference_lookup
            else np.nan
        )
        hq_is_missing = not np.isfinite(hq_estimate)
        all_is_missing = not np.isfinite(all_estimate)
        return (
            hq_is_missing,
            -float(hq_estimate) if not hq_is_missing else 0.0,
            all_is_missing,
            -float(all_estimate) if not all_is_missing else 0.0,
            str(amplicon),
        )

    amplicons = results["amplicon"].drop_duplicates().astype(str).tolist()
    return sorted(amplicons, key=sort_key)


def write_editing_rate_fixed_cell_depth_stability_plot(
    results: pd.DataFrame,
    output_root: str,
) -> List[Dict[str, str]]:
    """Write an all-amplicon exploratory stability plot at fixed cell counts."""
    plot_suffix = ".15_EditingRateFixedCellDepthStability"
    _remove_plot_artifacts(output_root, (plot_suffix,))
    if results.empty:
        return []

    all_amplicon_order = _fixed_cell_depth_stability_amplicon_order(results)
    numeric_columns = [
        "requested_sample_n_cells",
        "median_deviation_pp",
        "deviation_interval_lower_pp",
        "deviation_interval_upper_pp",
    ]
    plotted = results.copy()
    for column in numeric_columns:
        plotted[column] = pd.to_numeric(plotted[column], errors="coerce")
    plotted = plotted.loc[plotted["status"] == "ok"].dropna(
        subset=[
            "requested_sample_n_cells",
            "median_deviation_pp",
            "deviation_interval_lower_pp",
            "deviation_interval_upper_pp",
        ]
    )

    figure = _create_depth_stability_figure(
        plotted,
        all_amplicon_order,
        median_column="median_deviation_pp",
        lower_column="deviation_interval_lower_pp",
        upper_column="deviation_interval_upper_pp",
        title="Editing-Rate Stability Across Fixed Cell Counts",
        subtitle=(
            "Exploratory analysis; all amplicons are shown regardless of "
            "editing-rate significance"
        ),
        y_label="Deviation from full-cohort editing rate (percentage points)",
        x_column="requested_sample_n_cells",
        x_label="Eligible cells sampled",
        x_limits=(0.0, 1050.0),
        x_tick_values=FIXED_CELL_COUNT_DEPTHS,
        annotate_missing=True,
    )
    plot_root = output_root + plot_suffix
    _save_depth_stability_figure(figure, plot_root)
    return [
        {
            "plot_name": plot_root,
            "plot_title": "Editing-rate stability at fixed cell counts",
            "plot_label": (
                "Exploratory fixed-count downsampling bands for all analyzable "
                "and configured analysis-group cells at 100, 200, 400, 800, "
                "and 1,000 cells. Every amplicon is shown regardless of "
                "editing-rate significance; unavailable cohort sizes are not "
                "resampled or capped."
            ),
        }
    ]
