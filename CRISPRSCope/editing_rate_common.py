"""Shared editing-rate configuration, schemas, and numeric helpers."""

from dataclasses import dataclass
from typing import Sequence, Tuple

import numpy as np

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
    ".11_EditingRateCoverageAdjustedEffects",
    ".12_EditingRateUnconditionalPermutation",
    ".13_EditingRateObservedCenteredPermutationSwarm",
)

DEPTH_STABILITY_PLOT_SUFFIXES = (
    ".14_EditingRateDepthStability",
)

# Remove plots produced by earlier feature-branch versions so reruns cannot
# leave stale figures linked or mistaken for current output.
RETIRED_EDITING_RATE_PLOT_SUFFIXES = (
    ".11_EditingRateQualityDelta",
    ".12_EditingRateDepthStability",
    ".13_EditingRateRelativeDepthStability",
    ".14_EditingRateCoverageAdjustedEffects",
    ".15_EditingRateFixedCellDepthStability",
    ".16_EditingRateUnconditionalPermutation",
    ".17_EditingRateObservedCenteredPermutationSwarm",
    ".18_EditingRateDepthStability",
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

    enabled: bool = False
    iterations: int = 1_000
    percentages: Tuple[float, ...] = (10.0, 25.0, 50.0, 75.0, 90.0)
    confidence_level: float = 0.95
    seed: int = 42


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
