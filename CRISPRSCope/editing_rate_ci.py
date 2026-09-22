"""Compatibility façade for the modular editing-rate implementation.

Existing callers may continue importing calculations, configuration, and
test-used helpers from :mod:`CRISPRSCope.editing_rate_ci`.
"""

from .editing_rate_common import (
    CI_PLOT_SUFFIXES,
    DEPTH_STABILITY_COLUMNS,
    DEPTH_STABILITY_PLOT_SUFFIXES,
    RESULT_COLUMNS,
    RETIRED_EDITING_RATE_PLOT_SUFFIXES,
    UNCONDITIONAL_PERMUTATION_SIMULATION_COLUMNS,
    UNCONDITIONAL_PERMUTATION_SUMMARY_COLUMNS,
    EditingRateCIConfig,
    EditingRateDepthStabilityConfig,
    _percentile_interval,
    _point_estimate,
)
from .editing_rate_resampling import (
    _benjamini_hochberg,
    _bootstrap_means,
    _bootstrap_one_amplicon,
    _bootstrap_stratified_delta,
    _compute_editing_rate_resampling_analyses,
    _coverage_adjusted_resampling,
    _coverage_bin_ids,
    _permutation_hq_minus_all,
    _sample_empirical_means,
    _sample_without_replacement_sums,
    _two_sided_permutation_p_value,
    compute_editing_rate_confidence_intervals,
    compute_editing_rate_resampling_analyses,
)
from .editing_rate_depth_stability import (
    _depth_stability_one_amplicon,
    _empty_depth_stability_row,
    _nested_subsample_means,
    compute_editing_rate_depth_stability,
)
from . import editing_rate_plots as _plots
from .editing_rate_plots import (
    _create_depth_stability_figure,
    _depth_stability_amplicon_order,
    _depth_stability_layout,
    _finite_interval_rows,
    _interval_axis_limits,
    _observed_centered_permutation_differences,
    _observed_centered_swarm_coordinates,
    _remove_plot_artifacts,
    _row_plot_height,
    _save_depth_stability_figure,
    _significant_plot_rows,
    _write_coverage_adjusted_effect_plot,
    remove_editing_rate_plot_artifacts,
    write_editing_rate_ci_plots,
    write_editing_rate_observed_centered_permutation_swarm_plot,
    write_editing_rate_unconditional_permutation_plot,
)


def write_editing_rate_depth_stability_plot(*args, **kwargs):
    """Dispatch while preserving legacy façade monkeypatch behavior."""
    _plots._create_depth_stability_figure = _create_depth_stability_figure
    return _plots.write_editing_rate_depth_stability_plot(*args, **kwargs)


__all__ = [
    "CI_PLOT_SUFFIXES",
    "DEPTH_STABILITY_COLUMNS",
    "DEPTH_STABILITY_PLOT_SUFFIXES",
    "RESULT_COLUMNS",
    "RETIRED_EDITING_RATE_PLOT_SUFFIXES",
    "UNCONDITIONAL_PERMUTATION_SIMULATION_COLUMNS",
    "UNCONDITIONAL_PERMUTATION_SUMMARY_COLUMNS",
    "EditingRateCIConfig",
    "EditingRateDepthStabilityConfig",
    "_benjamini_hochberg",
    "_bootstrap_means",
    "_bootstrap_one_amplicon",
    "_bootstrap_stratified_delta",
    "_compute_editing_rate_resampling_analyses",
    "_coverage_adjusted_resampling",
    "_coverage_bin_ids",
    "_create_depth_stability_figure",
    "_depth_stability_amplicon_order",
    "_depth_stability_layout",
    "_depth_stability_one_amplicon",
    "_empty_depth_stability_row",
    "_finite_interval_rows",
    "_interval_axis_limits",
    "_nested_subsample_means",
    "_observed_centered_permutation_differences",
    "_observed_centered_swarm_coordinates",
    "_percentile_interval",
    "_permutation_hq_minus_all",
    "_point_estimate",
    "_remove_plot_artifacts",
    "_row_plot_height",
    "_sample_empirical_means",
    "_sample_without_replacement_sums",
    "_save_depth_stability_figure",
    "_significant_plot_rows",
    "_two_sided_permutation_p_value",
    "_write_coverage_adjusted_effect_plot",
    "compute_editing_rate_confidence_intervals",
    "compute_editing_rate_depth_stability",
    "compute_editing_rate_resampling_analyses",
    "remove_editing_rate_plot_artifacts",
    "write_editing_rate_ci_plots",
    "write_editing_rate_depth_stability_plot",
    "write_editing_rate_observed_centered_permutation_swarm_plot",
    "write_editing_rate_unconditional_permutation_plot",
]
