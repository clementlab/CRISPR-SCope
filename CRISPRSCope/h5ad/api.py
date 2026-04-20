"""High-level API for building h5ad output from CRISPRSCope pipeline artifacts."""

from pathlib import Path
from typing import Any, Dict, Optional
import logging
import shutil

from anndata import AnnData

from .builder import CRISPRSCopeAnnDataBuilder
from .loaders import (
    load_amplicons,
    load_crispresso_alleles,
    load_editing_summary,
    load_quality_scores,
    load_settings,
)
from .paths import resolve_h5ad_input_paths
from .utils import resolve_settings_relative_path

logger = logging.getLogger(__name__)


def get_default_h5ad_export_config() -> Dict[str, Any]:
    """Return the built-in h5ad export defaults."""
    return {
        "analysis_parameters": {
            "zygosity": {
                "wt_max_mod_pct": 20.0,
                "het_max_mod_pct": 80.0,
                "hom_min_mod_pct": 80.0,
                "compound_het_min_allele2_pct": 20.0,
            }
        }
    }


def build_h5ad_from_output_root(
    output_root: str,
    output_path: Optional[str] = None,
    config: Optional[Dict[str, Any]] = None,
    settings_path: Optional[str] = None,
    crispresso_dir: Optional[str] = None,
    n_processes: Optional[int] = None,
) -> AnnData:
    """Build and optionally save an AnnData object from CRISPRSCope pipeline outputs."""
    logger.info("Starting h5ad export for output_root=%s", output_root)
    export_config = config if config is not None else get_default_h5ad_export_config()
    paths = resolve_h5ad_input_paths(
        output_root=output_root,
        settings_path=settings_path,
        crispresso_dir=crispresso_dir,
    )

    allele_paths = []
    try:
        settings = load_settings(paths.settings_path)
        amplicons_value = settings.get("amplicons")
        if not amplicons_value:
            raise KeyError("The settings file does not contain an 'amplicons' entry.")

        amplicons_path = resolve_settings_relative_path(amplicons_value, paths.settings_path)
        amplicons_df = load_amplicons(amplicons_path)
        summary_df = load_editing_summary(paths.editing_summary_path)
        scores_df = load_quality_scores(paths.quality_scores_path)
        allele_paths = load_crispresso_alleles(paths.crispresso_filtered_dir, n_processes=n_processes)

        builder = CRISPRSCopeAnnDataBuilder(
            config=export_config,
            settings=settings,
            amplicons=amplicons_df,
            editing_summary=summary_df,
            quality_scores=scores_df,
            allele_parquet_paths=allele_paths,
        )
        adata = builder.build()

        if output_path:
            output_file = Path(output_path)
            output_file.parent.mkdir(parents=True, exist_ok=True)
            adata.write_h5ad(output_file, compression="gzip")
            logger.info("Saved h5ad output to %s", output_file)

        return adata
    finally:
        if allele_paths:
            temp_dir_to_remove = allele_paths[0].parent
            try:
                shutil.rmtree(temp_dir_to_remove)
                logger.info("Removed temporary allele directory %s", temp_dir_to_remove)
            except Exception as exc:
                logger.warning("Could not remove temporary allele directory %s: %s", temp_dir_to_remove, exc)
