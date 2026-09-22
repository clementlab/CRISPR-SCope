"""Plotting and stale-artifact cleanup for editing-rate analyses."""

import logging
import math
import os
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from CRISPRSCope.output_artifacts import OutputContext

from .editing_rate_common import (
    CI_PLOT_SUFFIXES,
    DEPTH_STABILITY_PLOT_SUFFIXES,
    RETIRED_EDITING_RATE_PLOT_SUFFIXES,
)

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
    outputs = OutputContext(output_root)
    for suffix in suffixes:
        registered_key = outputs.key_for_suffix(suffix)
        if registered_key is not None:
            outputs.remove_optional((registered_key,))
            continue
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
        CI_PLOT_SUFFIXES
        + DEPTH_STABILITY_PLOT_SUFFIXES
        + RETIRED_EDITING_RATE_PLOT_SUFFIXES,
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

    plotted["_coverage_adjusted_distance_from_zero"] = np.abs(
        plotted["coverage_adjusted_in_group_minus_out_group_pct"]
    )
    plotted = plotted.sort_values(
        ["_coverage_adjusted_distance_from_zero", "amplicon"],
        ascending=[False, True],
        na_position="last",
        kind="mergesort",
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
    ax.invert_yaxis()
    ax.set_xlabel("InGroup minus OutGroup (percentage points)")
    ax.set_ylabel("Amplicon")
    fig.suptitle("Raw and Coverage-Adjusted Editing-Rate Effects")
    if legend_handles:
        ax.legend(
            [item[0] for item in legend_handles],
            [item[1] for item in legend_handles],
            title="Open = not significant; filled = corresponding BH-adjusted p-value ≤ 0.05",
            loc="lower center",
            bbox_to_anchor=(0.5, 1.01),
            borderaxespad=0.0,
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
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.90))
    plot_root = OutputContext(output_root).plot_root(
        "editing_rate_coverage_adjusted_effects_plot"
    )
    fig.savefig(plot_root + ".pdf", bbox_inches="tight")
    fig.savefig(plot_root + ".png", bbox_inches="tight")
    plt.close(fig)
    return [
        {
            "artifact_key": "editing_rate_coverage_adjusted_effects_plot",
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
    suffix = ".12_EditingRateUnconditionalPermutation"
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
    plot_root = OutputContext(output_root).plot_root(
        "editing_rate_unconditional_permutation_plot"
    )
    fig.savefig(plot_root + ".pdf", bbox_inches="tight")
    fig.savefig(plot_root + ".png", bbox_inches="tight")
    plt.close(fig)
    return [
        {
            "artifact_key": "editing_rate_unconditional_permutation_plot",
            "plot_name": plot_root,
            "plot_title": "Unconditional editing-rate permutation distribution",
            "plot_label": (
                "For each amplicon, boxplots show configured-group-sized subsets "
                "drawn without replacement from all eligible cells. Orange diamonds "
                "show the observed configured analysis-group estimates."
            ),
        }
    ]


def _observed_centered_permutation_differences(
    summaries: pd.DataFrame,
    simulations: pd.DataFrame,
) -> pd.DataFrame:
    """Return valid unconditional draws centered on each observed group mean."""
    required_summary_columns = {
        "amplicon",
        "selected_group_estimate_pct",
        "valid_permutations",
        "seed",
    }
    required_simulation_columns = {"amplicon", "permuted_selected_estimate_pct"}
    if (
        not required_summary_columns.issubset(summaries.columns)
        or not required_simulation_columns.issubset(simulations.columns)
    ):
        return pd.DataFrame(
            columns=(
                "amplicon",
                "observed_in_group_estimate_pct",
                "simulated_minus_observed_pct",
                "seed",
            )
        )

    observed = summaries.loc[:, [
        "amplicon",
        "selected_group_estimate_pct",
        "valid_permutations",
        "seed",
    ]].copy()
    observed["selected_group_estimate_pct"] = pd.to_numeric(
        observed["selected_group_estimate_pct"], errors="coerce"
    )
    observed["valid_permutations"] = pd.to_numeric(
        observed["valid_permutations"], errors="coerce"
    )
    observed = observed.loc[
        observed["selected_group_estimate_pct"].notna()
        & observed["valid_permutations"].gt(0)
    ].drop_duplicates("amplicon")
    if observed.empty:
        return pd.DataFrame()

    draws = simulations.loc[:, ["amplicon", "permuted_selected_estimate_pct"]].copy()
    draws["permuted_selected_estimate_pct"] = pd.to_numeric(
        draws["permuted_selected_estimate_pct"], errors="coerce"
    )
    centered = draws.merge(observed, on="amplicon", how="inner", validate="many_to_one")
    centered = centered.loc[centered["permuted_selected_estimate_pct"].notna()].copy()
    centered["simulated_minus_observed_pct"] = (
        centered["permuted_selected_estimate_pct"]
        - centered["selected_group_estimate_pct"]
    )
    return centered.rename(
        columns={"selected_group_estimate_pct": "observed_in_group_estimate_pct"}
    )[[
        "amplicon",
        "observed_in_group_estimate_pct",
        "simulated_minus_observed_pct",
        "seed",
    ]]


def _observed_centered_swarm_coordinates(
    centered: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Order draws by mean relative effect and assign reproducible vertical jitter."""
    observed = centered.loc[:, ["amplicon", "observed_in_group_estimate_pct", "seed"]]
    observed = observed.drop_duplicates("amplicon")
    mean_relative_effect = (
        centered.groupby("amplicon", as_index=False)["simulated_minus_observed_pct"]
        .mean()
        .rename(columns={"simulated_minus_observed_pct": "mean_relative_effect_pct"})
    )
    observed = observed.merge(mean_relative_effect, on="amplicon", validate="one_to_one")
    observed = observed.sort_values(
        ["mean_relative_effect_pct", "amplicon"],
        ascending=[False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    order = {amplicon: index for index, amplicon in enumerate(observed["amplicon"])}
    coordinates = centered.copy().reset_index(drop=True)
    coordinates["row_index"] = coordinates["amplicon"].map(order)
    jitter = np.empty(len(coordinates), dtype=float)
    for amplicon, group_indices in coordinates.groupby("amplicon", sort=False).groups.items():
        row = observed.loc[observed["amplicon"] == amplicon].iloc[0]
        seed = int(pd.to_numeric(row["seed"], errors="coerce"))
        rng = np.random.default_rng(np.random.SeedSequence([seed, order[amplicon]]))
        jitter[list(group_indices)] = rng.uniform(-0.16, 0.16, size=len(group_indices))
    coordinates["swarm_y"] = coordinates["row_index"].to_numpy(dtype=float) + jitter
    return coordinates, observed


def write_editing_rate_observed_centered_permutation_swarm_plot(
    summaries: pd.DataFrame,
    simulations: pd.DataFrame,
    output_root: str,
) -> List[Dict[str, str]]:
    """Plot unconditional draws as differences from observed InGroup means."""
    suffix = ".13_EditingRateObservedCenteredPermutationSwarm"
    _remove_plot_artifacts(output_root, (suffix,))
    centered = _observed_centered_permutation_differences(summaries, simulations)
    if centered.empty:
        return []
    centered, observed = _observed_centered_swarm_coordinates(centered)

    fig, ax = plt.subplots(figsize=(12.0, _row_plot_height(len(observed))))
    ax.scatter(
        centered["simulated_minus_observed_pct"],
        centered["swarm_y"],
        s=7,
        color="#4C78A8",
        alpha=0.25,
        linewidths=0,
        rasterized=True,
    )
    ax.axvline(
        0,
        color="black",
        linestyle="--",
        linewidth=1,
        label="Observed InGroup editing rate",
    )
    ax.set_yticks(np.arange(len(observed), dtype=float))
    ax.set_yticklabels(observed["amplicon"])
    ax.set_xlabel("Simulated mean relative to observed InGroup mean (percentage points)")
    ax.set_ylabel("Amplicon")
    ax.set_title("Observed-Centered Unconditional Permutation Swarm")
    limits = _interval_axis_limits(
        centered["simulated_minus_observed_pct"].to_numpy(dtype=float),
        centered["simulated_minus_observed_pct"].to_numpy(dtype=float),
        minimum_span=2.0,
        include_zero=True,
    )
    if limits is not None:
        ax.set_xlim(*limits)
    ax.grid(axis="x", alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    plot_root = OutputContext(output_root).plot_root(
        "editing_rate_observed_centered_permutation_swarm_plot"
    )
    fig.savefig(plot_root + ".pdf", bbox_inches="tight")
    fig.savefig(plot_root + ".png", bbox_inches="tight")
    plt.close(fig)
    return [
        {
            "artifact_key": "editing_rate_observed_centered_permutation_swarm_plot",
            "plot_name": plot_root,
            "plot_title": "Observed-centered unconditional permutation swarm",
            "plot_label": (
                "Each point is an InGroup-sized subset drawn without replacement "
                "from eligible cells, shown relative to the observed InGroup "
                "editing rate for that amplicon; zero marks the observed rate."
            ),
        }
    ]


def write_editing_rate_ci_plots(results: pd.DataFrame, output_root: str) -> List[Dict[str, str]]:
    """Plot controlled-significant intervals and all adjusted-effect comparisons."""
    plot_metadata: List[Dict[str, str]] = []
    _remove_plot_artifacts(
        output_root,
        CI_PLOT_SUFFIXES + RETIRED_EDITING_RATE_PLOT_SUFFIXES,
    )
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
        plot_root = OutputContext(output_root).plot_root(
            "editing_rate_confidence_intervals_plot"
        )
        fig.savefig(plot_root + ".pdf", bbox_inches="tight")
        fig.savefig(plot_root + ".png", bbox_inches="tight")
        plt.close(fig)
        plot_metadata.append(
            {
                "artifact_key": "editing_rate_confidence_intervals_plot",
                "plot_name": plot_root,
                "plot_title": "Amplicon editing-rate confidence intervals",
                "plot_label": (
                    "Pointwise bootstrap confidence intervals for all analyzable "
                    "and configured analysis-group cells among amplicons with a "
                    "coverage-adjusted BH p-value at or below 0.05."
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
    """Write the absolute stability plot with optional significance filtering."""
    _remove_plot_artifacts(
        output_root,
        DEPTH_STABILITY_PLOT_SUFFIXES + RETIRED_EDITING_RATE_PLOT_SUFFIXES,
    )
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
    absolute_plot_root = OutputContext(output_root).plot_root(
        "editing_rate_depth_stability_plot"
    )
    _save_depth_stability_figure(absolute_figure, absolute_plot_root)
    significance_clause = (
        " among amplicons with a coverage-adjusted BH p-value at or below 0.05"
        if significance_filter_applied
        else ""
    )
    return [
        {
            "artifact_key": "editing_rate_depth_stability_plot",
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
