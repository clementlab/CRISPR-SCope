"""Path resolution helpers for h5ad export."""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass(frozen=True)
class H5ADInputPaths:
    """Resolved CRISPRSCope output paths needed for h5ad export."""

    output_root: Path
    settings_path: Path
    editing_summary_path: Path
    quality_scores_path: Path
    crispresso_filtered_dir: Path


def default_h5ad_output_path(output_root: str) -> Path:
    """Return the default .h5ad path for a pipeline output root."""
    return Path(f"{output_root}.h5ad")


def resolve_h5ad_input_paths(
    output_root: str,
    settings_path: Optional[str] = None,
    crispresso_dir: Optional[str] = None,
) -> H5ADInputPaths:
    """Resolve the pipeline artifacts needed to build an h5ad export."""
    root = Path(output_root)

    if settings_path is not None:
        resolved_settings = Path(settings_path)
    elif root.is_file():
        resolved_settings = root
    else:
        raise FileNotFoundError(
            "Unable to infer the original settings file from output_root. "
            "Provide settings_path explicitly when output_root differs from the settings file path."
        )

    return H5ADInputPaths(
        output_root=root,
        settings_path=resolved_settings,
        editing_summary_path=Path(f"{output_root}.filteredEditingSummary.txt"),
        quality_scores_path=Path(f"{output_root}.amplicon_score.txt"),
        crispresso_filtered_dir=Path(crispresso_dir) if crispresso_dir else Path(f"{output_root}.crispresso.filtered"),
    )
