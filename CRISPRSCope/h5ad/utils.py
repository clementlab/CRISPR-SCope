"""Utility helpers for h5ad export."""

from pathlib import Path
from typing import List
import logging

import pandas as pd

logger = logging.getLogger(__name__)


def extract_amplicon_names(editing_df: pd.DataFrame) -> List[str]:
    """Extract unique amplicon names from editing-summary columns."""
    amplicons = set()
    for col in editing_df.columns:
        if col.startswith(("totCount.", "modPct.")):
            amplicons.add(col.split(".", 1)[1])

    amplicon_list = sorted(amplicons)
    logger.info("Extracted %d unique amplicon names.", len(amplicon_list))
    return amplicon_list


def resolve_settings_relative_path(path_value: str, settings_path: Path) -> Path:
    """Resolve a settings-derived path relative to the settings file directory."""
    candidate = Path(path_value)
    if candidate.is_absolute():
        return candidate
    return (settings_path.parent / candidate).resolve()
