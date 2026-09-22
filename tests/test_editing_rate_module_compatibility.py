"""Ownership and legacy-import guarantees for modular editing-rate code."""

from CRISPRSCope import editing_rate_ci, pipeline
from CRISPRSCope import (
    editing_rate_common,
    editing_rate_depth_stability,
    editing_rate_outputs,
    editing_rate_plots,
    editing_rate_resampling,
)


def test_editing_rate_ci_facade_reexports_owner_symbols():
    assert editing_rate_ci.EditingRateCIConfig is editing_rate_common.EditingRateCIConfig
    assert editing_rate_ci.RESULT_COLUMNS is editing_rate_common.RESULT_COLUMNS
    assert editing_rate_ci._bootstrap_means is editing_rate_resampling._bootstrap_means
    assert (
        editing_rate_ci.compute_editing_rate_resampling_analyses
        is editing_rate_resampling.compute_editing_rate_resampling_analyses
    )
    assert (
        editing_rate_ci.compute_editing_rate_depth_stability
        is editing_rate_depth_stability.compute_editing_rate_depth_stability
    )
    assert (
        editing_rate_ci.write_editing_rate_ci_plots
        is editing_rate_plots.write_editing_rate_ci_plots
    )
    assert (
        editing_rate_ci.remove_editing_rate_plot_artifacts
        is editing_rate_plots.remove_editing_rate_plot_artifacts
    )


def test_pipeline_reexports_relocated_editing_rate_output_writers():
    assert pipeline.write_editing_rate_ci_output is editing_rate_outputs.write_editing_rate_ci_output
    assert (
        pipeline.write_editing_rate_depth_stability_output
        is editing_rate_outputs.write_editing_rate_depth_stability_output
    )
