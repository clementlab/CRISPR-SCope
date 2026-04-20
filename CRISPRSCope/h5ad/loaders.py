"""Data loading helpers for h5ad export."""

from pathlib import Path
from typing import Any, Dict, List, Optional
import gzip
import logging
import tempfile
import traceback
from collections import Counter
from multiprocessing import Pool, cpu_count

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

logger = logging.getLogger(__name__)


def _parse_and_write_parquet(input_fastq_path: Path, output_parquet_path: Path) -> Optional[str]:
    """Parse a single CRISPResso_output.fastq.gz file and persist allele counts to Parquet."""
    amplicon_name = input_fastq_path.parent.name.replace("CRISPResso_on_", "", 1)
    allele_counts = Counter()

    try:
        with gzip.open(input_fastq_path, "rt") as handle:
            for header in handle:
                sequence = next(handle).strip()
                next(handle)
                next(handle)
                cell_barcode = header.strip().split(":")[-2]
                allele_counts[(cell_barcode, sequence)] += 1
    except StopIteration:
        logger.warning("File %s appears to be truncated or empty.", input_fastq_path)
        return None
    except Exception as exc:
        logger.error("Error parsing file %s: %s\n%s", input_fastq_path, exc, traceback.format_exc())
        return None

    if not allele_counts:
        return None

    results = [
        {
            "cell_barcode": cell_barcode,
            "amplicon_name": amplicon_name,
            "allele_sequence": sequence,
            "count": count,
        }
        for (cell_barcode, sequence), count in allele_counts.items()
    ]

    try:
        df = pd.DataFrame(results)
        df["cell_barcode"] = df["cell_barcode"].astype("category")
        df["amplicon_name"] = df["amplicon_name"].astype("category")
        table = pa.Table.from_pandas(df, preserve_index=False)
        pq.write_table(table, output_parquet_path, compression="snappy")
        return str(output_parquet_path)
    except Exception as exc:
        logger.error("Failed to write Parquet file %s: %s", output_parquet_path, exc)
        return None


def load_settings(settings_path: Path) -> Dict[str, Any]:
    """Parse a CRISPRSCope settings file."""
    if not settings_path.is_file():
        logger.error("Settings file not found at: %s", settings_path)
        raise FileNotFoundError(f"Settings file not found: {settings_path}")

    settings: Dict[str, Any] = {}
    with open(settings_path, "r") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            if len(parts) == 2:
                key, value = parts
                settings[key] = value
            else:
                logger.warning("Could not parse line in settings file: %r", line)
    logger.info("Loaded %d parameters from %s.", len(settings), settings_path)
    return settings


def load_amplicons(amplicons_path: Path) -> pd.DataFrame:
    """Load amplicon metadata."""
    if not amplicons_path.is_file():
        logger.error("Amplicons file not found at: %s", amplicons_path)
        raise FileNotFoundError(f"Amplicons file not found: {amplicons_path}")

    df = pd.read_csv(
        amplicons_path,
        sep="\t",
        header=None,
        index_col=0,
        names=["sequence", "guide"],
    )
    df.index.name = "amplicon_name"
    logger.info("Loaded %d amplicons from %s.", len(df), amplicons_path)
    return df


def load_editing_summary(summary_path: Path) -> pd.DataFrame:
    """Load the filtered editing summary."""
    if not summary_path.is_file():
        logger.error("Editing summary file not found at: %s", summary_path)
        raise FileNotFoundError(f"Editing summary file not found: {summary_path}")

    df = pd.read_csv(summary_path, sep="\t", index_col=0)
    df.index.name = "cell_barcode"
    logger.info("Loaded editing data for %d cells and %d metrics.", df.shape[0], df.shape[1])
    return df


def load_quality_scores(scores_path: Path) -> pd.DataFrame:
    """Load cell-level quality annotations."""
    if not scores_path.is_file():
        logger.error("Amplicon scores file not found at: %s", scores_path)
        raise FileNotFoundError(f"Amplicon scores file not found: {scores_path}")

    df = pd.read_csv(scores_path, sep="\t", index_col=0)
    df.index.name = "cell_barcode"
    logger.info("Loaded quality scores for %d cells.", len(df))
    return df


def load_crispresso_alleles(crispresso_dir: Path, n_processes: Optional[int] = None) -> List[Path]:
    """Parse filtered CRISPResso FASTQ outputs into persistent intermediate Parquet files."""
    if not crispresso_dir.is_dir():
        logger.error("CRISPResso output directory not found at: %s", crispresso_dir)
        raise FileNotFoundError(f"Directory not found: {crispresso_dir}")

    input_files = list(crispresso_dir.glob("CRISPResso_on_*/CRISPResso_output.fastq.gz"))
    if not input_files:
        logger.warning("No CRISPResso FASTQ files found in %s.", crispresso_dir)
        return []

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        job_args = [
            (input_path, temp_path / f"{input_path.parent.name}.parquet")
            for input_path in input_files
        ]

        worker_count = n_processes or cpu_count()
        logger.info("Parsing %d CRISPResso FASTQs with %d processes.", len(job_args), worker_count)

        with Pool(processes=worker_count) as pool:
            results = pool.starmap(_parse_and_write_parquet, job_args)
            temp_parquet_paths = [Path(path) for path in results if path is not None]

        if not temp_parquet_paths:
            logger.warning("Allele parsing yielded no intermediate Parquet files.")
            return []

        persistent_temp_dir = Path(tempfile.mkdtemp(prefix="CRISPRSCope_h5ad_alleles_"))
        final_paths = []
        for temp_file in temp_parquet_paths:
            final_path = persistent_temp_dir / temp_file.name
            temp_file.rename(final_path)
            final_paths.append(final_path)

    logger.info("Prepared %d allele Parquet files in %s.", len(final_paths), persistent_temp_dir)
    return final_paths
