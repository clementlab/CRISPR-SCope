"""Validate guide placement before the expensive read-processing stages."""
from __future__ import annotations

import logging
import os
from dataclasses import dataclass

import numpy as np
from CRISPResso2 import CRISPResso2Align, CRISPRessoShared


FLEXIGUIDE_HOMOLOGY_DEFAULT = 80


@dataclass(frozen=True)
class GuideResolution:
    mode: str
    exact: tuple[str, ...]
    flexible: tuple[str, ...]
    omitted: tuple[str, ...]
    homology: int


def _matches_exact(guide: str, references: list[str]) -> bool:
    return any(
        guide in reference or CRISPRessoShared.reverse_complement(guide) in reference
        for reference in references
    )


def _matches_flexible(guide: str, references: list[str], homology: int, matrix) -> bool:
    """Mirror the installed CRISPResso2 flexiguide search and strict score cutoff."""
    for reference in references:
        for oriented_ref in (reference, CRISPRessoShared.reverse_complement(reference)):
            incentive = np.zeros(len(oriented_ref) + 1, dtype=int)
            aligned_guide, _, _ = CRISPResso2Align.global_align(
                guide, oriented_ref, matrix=matrix, gap_incentive=incentive,
                gap_open=-20, gap_extend=-2,
            )
            candidate = aligned_guide.strip("-")
            if abs(len(candidate) - len(guide)) >= 2:
                continue
            position = aligned_guide.find(candidate)
            locus = oriented_ref[position:position + len(candidate)]
            incentive = np.zeros(len(locus) + 1, dtype=int)
            _, _, score = CRISPResso2Align.global_align(
                guide, locus, matrix=matrix, gap_incentive=incentive,
                gap_open=-20, gap_extend=-2,
            )
            if score > homology:
                return True
    return False


def validate_guides(amplicon_file: str, alternate_alleles_file: str = "",
                    homology: int = FLEXIGUIDE_HOMOLOGY_DEFAULT) -> dict[str, GuideResolution]:
    """Resolve all guides, reporting exact errors together before processing reads."""
    if not isinstance(homology, int) or not 0 <= homology <= 100:
        raise ValueError("flexiguide_homology must be an integer between 0 and 100")

    alternate_refs = {}
    if alternate_alleles_file:
        with open(alternate_alleles_file) as handle:
            next(handle, None)  # Header; same layout used by the amplicon splitter.
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 3:
                    alternate_refs[fields[0]] = fields[2]

    rows = []
    with open(amplicon_file) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"Amplicon file line {line_number} requires a name and sequence")
            name = fields[0].strip()
            reference_text = alternate_refs.get(name, fields[1])
            references = [sequence.strip().upper() for sequence in reference_text.split(",")]
            raw_guide = fields[2].strip() if len(fields) > 2 else ""
            mode = fields[4].strip().lower() if len(fields) > 4 else ""
            mode = mode or "exact"
            if mode not in ("exact", "flexible"):
                raise ValueError(
                    f"Amplicon {name} (line {line_number}): guide_match_mode must be exact or flexible, got {mode!r}"
                )
            guides = [] if raw_guide.lower() in ("", "na", "none") else [
                guide.strip().upper() for guide in raw_guide.split(",")
            ]
            rows.append((line_number, name, references, guides, mode))

    # Loading the matrix is inexpensive, but unnecessary when there are no flexible mismatches.
    matrix = None
    resolutions = {}
    errors = []
    for line_number, name, references, guides, mode in rows:
        exact, flexible, omitted = [], [], []
        for guide in guides:
            if _matches_exact(guide, references):
                exact.append(guide)
            elif mode == "exact":
                errors.append(
                    f"  - Amplicon {name} (line {line_number}): guide RNA {guide} does not perfectly "
                    "match any reference sequence on either strand. Check the sequence or set "
                    "guide_match_mode to flexible."
                )
            else:
                if matrix is None:
                    matrix_path = os.path.join(os.path.dirname(CRISPRessoShared.__file__), "EDNAFULL")
                    matrix = CRISPResso2Align.read_matrix(matrix_path)
                if _matches_flexible(guide, references, homology, matrix):
                    flexible.append(guide)
                else:
                    omitted.append(guide)
        for guide in omitted:
            logging.warning(
                "Amplicon %s (line %d): guide RNA %s did not exceed flexiguide_homology %d; "
                "omitting this guide%s.", name, line_number, guide, homology,
                " (CRISPResso2 will quantify the full amplicon)" if not exact and not flexible else "",
            )
        resolutions[name] = GuideResolution(mode, tuple(exact), tuple(flexible), tuple(omitted), homology)
    if errors:
        message = "Guide validation failed before read processing:\n" + "\n".join(errors)
        logging.error("%s", message)
        raise ValueError(message)
    return resolutions
