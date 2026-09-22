"""Finite-cohort editing-rate depth-stability analysis."""

import math
import multiprocessing as mp
from typing import Dict, Iterable, List, Sequence

import numpy as np
import pandas as pd

from .editing_rate_common import (
    DEPTH_STABILITY_COLUMNS,
    EditingRateDepthStabilityConfig,
    _percentile_interval,
    categorical_mod_pct_to_cell_edit_pct,
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
                "mod_values": categorical_mod_pct_to_cell_edit_pct(
                    pd.to_numeric(
                        ordered_summary[mod_column], errors="coerce"
                    ).to_numpy(dtype=float),
                    amplicon=amplicon,
                    barcodes=ordered_index,
                ),
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
    return pd.DataFrame(results).reindex(columns=DEPTH_STABILITY_COLUMNS)
