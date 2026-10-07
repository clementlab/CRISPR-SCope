"""Small, deterministic helpers for first-pass allele calling."""

from __future__ import annotations

import re
from dataclasses import dataclass
from statistics import median
from typing import Dict, List, Optional, Set, Union


@dataclass(frozen=True)
class AlleleSupport:
    mode: str
    value: Union[int, float]

    def allows(self, count: int, total: int) -> bool:
        if self.mode == "count":
            return count >= self.value
        return total > 0 and count / total >= self.value

    def cache_value(self) -> str:
        return f"{self.mode}:{self.value}"


def parse_allele_support(value: object) -> AlleleSupport:
    """An integer is a read count; decimal syntax is a fraction."""
    if isinstance(value, bool):
        raise ValueError("min_allele_support must be an integer count or decimal fraction")
    if isinstance(value, int) or isinstance(value, str) and re.fullmatch(r"\d+", value.strip()):
        count = int(value)
        return AlleleSupport("count", count)
    if isinstance(value, float) or isinstance(value, str) and re.fullmatch(r"(?:\d+\.\d*|\.\d+)", value.strip()):
        fraction = float(value)
        if not 0 <= fraction <= 1:
            raise ValueError("decimal min_allele_support must be between 0.0 and 1.0")
        return AlleleSupport("fraction", fraction)
    raise ValueError("min_allele_support must be an integer count (2) or decimal fraction (0.1)")


def annotation_fields(annotation: str) -> Optional[Dict[str, str]]:
    fields = {}
    for name in ("ALN", "DEL", "INS", "SUB", "ALN_REF", "ALN_SEQ"):
        match = re.search(r"(?:^|\s)" + name + r"=([^\s]*)", annotation)
        if match:
            fields[name] = match.group(1)
    return fields if fields.get("ALN") and "DEL" in fields and "INS" in fields else None


def allele_status(fields: Dict[str, str], ignore_substitutions: bool) -> Optional[str]:
    if not ignore_substitutions and "SUB" not in fields:
        return None
    status = f"DEL={fields['DEL']} INS={fields['INS']}"
    if not ignore_substitutions:
        status += f" SUB={fields['SUB']}"
    return status


def is_wildtype_key(key: str) -> bool:
    status = key.partition(":")[2]
    return status in ("DEL= INS=", "DEL= INS= SUB=")


def _edit_sites(fields: Dict[str, str], ignore_substitutions: bool) -> Optional[Dict[str, Set[int]]]:
    """Decode CRISPResso's quantification-window edit positions."""
    sites = {"DEL": set(), "INS": set(), "SUB": set()}
    for token in filter(None, fields["DEL"].split(";")):
        match = re.fullmatch(r"(\d+)\((\d+)\)", token)
        if match is None:
            return None
        start, size = (int(value) for value in match.groups())
        sites["DEL"].update(range(start, start + size))
    for token in filter(None, fields["INS"].split(";")):
        match = re.fullmatch(r"(\d+)\((\d+)\+([A-Za-z]+)\)", token)
        if match is None or int(match.group(2)) != len(match.group(3)):
            return None
        sites["INS"].add(int(match.group(1)))
    if not ignore_substitutions:
        for token in filter(None, fields.get("SUB", "").split(";")):
            if not token.isdigit():
                return None
            sites["SUB"].add(int(token))
    return sites


@dataclass
class LocalQuality:
    reference_qualities: Dict[int, int]
    event_sites: Set[int]
    event_score: Optional[int]

    def score_at_sites(self, sites: Set[int]) -> Optional[int]:
        values = [
            quality
            for site in sites
            for position, quality in self.reference_qualities.items()
            if abs(position - site) <= 2
        ]
        return min(values) if values else None


def local_quality(sequence: str, quality: str, fields: Dict[str, str], ignore_substitutions: bool) -> Optional[LocalQuality]:
    """Map FASTQ Phred scores through CRISPResso's gapped alignment."""
    aligned_ref = fields.get("ALN_REF")
    aligned_read = fields.get("ALN_SEQ")
    if (
        not aligned_ref or not aligned_read or len(aligned_ref) != len(aligned_read)
        or aligned_read.replace("-", "") != sequence or len(quality) != len(sequence)
    ):
        return None
    annotated_sites = _edit_sites(fields, ignore_substitutions)
    if annotated_sites is None:
        return None

    reference_qualities = {}
    events = []
    matched_sites = {"DEL": set(), "INS": set(), "SUB": set()}
    read_index = 0
    reference_index = 0
    for ref_base, read_base in zip(aligned_ref, aligned_read):
        ref_position = reference_index if ref_base != "-" else max(reference_index - 1, 0)
        phred = ord(quality[read_index]) - 33 if read_base != "-" else None
        if ref_base != "-" and phred is not None:
            reference_qualities[reference_index] = phred
        if ref_base == "-" and read_base != "-" and ref_position in annotated_sites["INS"]:
            events.append((ref_position, phred))
            matched_sites["INS"].add(ref_position)
        elif read_base == "-" and ref_base != "-" and ref_position in annotated_sites["DEL"]:
            events.append((ref_position, None))
            matched_sites["DEL"].add(ref_position)
        elif (
            ref_base != "-" and read_base != "-" and ref_base.upper() != read_base.upper()
            and ref_position in annotated_sites["SUB"]
        ):
            events.append((ref_position, phred))
            matched_sites["SUB"].add(ref_position)
        if read_base != "-":
            read_index += 1
        if ref_base != "-":
            reference_index += 1

    if any(matched_sites[kind] != annotated_sites[kind] for kind in matched_sites):
        return None

    sites = {position for position, _ in events}
    local_bases = [quality for _, quality in events if quality is not None]
    local_bases.extend(
        quality
        for site in sites
        for position, quality in reference_qualities.items()
        if abs(position - site) <= 2
    )
    return LocalQuality(reference_qualities, sites, min(local_bases) if local_bases else None)


def median_read_minimum(records: List[Optional[LocalQuality]], sites: Optional[Set[int]] = None) -> Optional[float]:
    """Return no score unless every supporting read can be scored."""
    if not records:
        return None
    scores = [record.score_at_sites(sites) if sites is not None else record.event_score for record in records if record is not None]
    if len(scores) != len(records) or any(score is None for score in scores):
        return None
    return float(median(scores))
