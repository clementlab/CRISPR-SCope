"""Versioned cache records for resumable CRISPRSCope processing stages.

The cache is intentionally small and declarative.  Pipeline stages declare the
inputs and parameters that determine their result, plus the files that prove a
completed result is still usable.  This module never infers dependencies from a
settings file and never hashes large sequencing artifacts during validation.
"""

from __future__ import annotations

import fcntl
import gzip
import hashlib
import json
import logging
import os
import re
import shutil
import socket
import struct
import subprocess
import tempfile
from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
from functools import lru_cache
from pathlib import Path
from typing import Callable, Iterable, Mapping, Sequence


CACHE_SCHEMA_VERSION = 1


class CacheMode(str, Enum):
    AUTO = "auto"
    REFRESH = "refresh"
    DISABLED = "disabled"


@dataclass(frozen=True)
class CacheConfig:
    mode: CacheMode = CacheMode.AUTO

    @classmethod
    def from_value(cls, value: object = "auto") -> "CacheConfig":
        normalized = str(value).strip().lower()
        try:
            return cls(mode=CacheMode(normalized))
        except ValueError as error:
            choices = ", ".join(mode.value for mode in CacheMode)
            raise ValueError(
                f"Invalid value for cache_mode: {value!r}. Expected one of: {choices}."
            ) from error


def _canonical_value(value: object) -> object:
    """Return a JSON-compatible value with deterministic mapping/set order."""
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        return {
            str(key): _canonical_value(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (set, frozenset)):
        normalized = [_canonical_value(item) for item in value]
        return sorted(
            normalized,
            key=lambda item: json.dumps(item, sort_keys=True, separators=(",", ":")),
        )
    if isinstance(value, (list, tuple)):
        return [_canonical_value(item) for item in value]
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    raise TypeError(f"Unsupported cache-key value: {type(value).__name__}")


def canonical_json(value: object) -> str:
    return json.dumps(
        _canonical_value(value), sort_keys=True, separators=(",", ":"), ensure_ascii=False
    )


def canonical_digest(value: object) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def _absolute_path(path: os.PathLike[str] | str) -> str:
    """Return one canonical spelling for a path, including symlink aliases."""
    return os.path.realpath(os.path.abspath(os.fspath(path)))


def small_file_fingerprint(path: os.PathLike[str] | str) -> dict[str, object]:
    """Fingerprint a small file by content."""
    absolute = _absolute_path(path)
    stat = os.stat(absolute)
    digest = hashlib.sha256()
    with open(absolute, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return {
        "path": absolute,
        "strategy": "sha256",
        "size": stat.st_size,
        "sha256": digest.hexdigest(),
    }


def large_file_fingerprint(path: os.PathLike[str] | str) -> dict[str, object]:
    """Fingerprint a large file without reading its payload."""
    absolute = _absolute_path(path)
    stat = os.stat(absolute)
    return {
        "path": absolute,
        "strategy": "stat",
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def gzip_content_fingerprint(path: os.PathLike[str] | str) -> dict[str, object]:
    """Fingerprint gzip content in constant time using its RFC 1952 trailer.

    Split FASTQs are regenerated as complete single-member gzip files.  Their
    CRC32 and uncompressed size therefore provide a stable content signature
    without rereading the sequencing payload.  Compressed size is retained as
    an additional guard against corruption or accidental substitution.
    """
    absolute = _absolute_path(path)
    stat = os.stat(absolute)
    if stat.st_size < 18:
        raise ValueError(f"Gzip file is too short to contain a valid trailer: {absolute}")
    with open(absolute, "rb") as handle:
        if handle.read(2) != b"\x1f\x8b":
            raise ValueError(f"File does not have a gzip signature: {absolute}")
        handle.seek(-8, os.SEEK_END)
        crc32, uncompressed_size = struct.unpack("<II", handle.read(8))
    return {
        "path": absolute,
        "strategy": "gzip_crc32",
        "size": stat.st_size,
        "crc32": f"{crc32:08x}",
        "uncompressed_size": uncompressed_size,
    }


def optional_file_fingerprint(
    path: os.PathLike[str] | str | None, *, strategy: str
) -> dict[str, object] | None:
    if not path:
        return None
    if strategy == "sha256":
        return small_file_fingerprint(path)
    if strategy == "stat":
        return large_file_fingerprint(path)
    raise ValueError(f"Unknown fingerprint strategy: {strategy}")


@lru_cache(maxsize=None)
def tool_identity(command: tuple[str, ...]) -> dict[str, str]:
    """Resolve and identify a tool once per Python process."""
    if not command:
        raise ValueError("Tool command must not be empty")
    executable = shutil.which(command[0])
    if executable is None:
        raise FileNotFoundError(f"Required executable not found: {command[0]}")
    output = subprocess.check_output(command, stderr=subprocess.STDOUT, text=True)
    if isinstance(output, bytes):
        output = output.decode("utf-8", errors="replace")
    normalized = "\n".join(line.rstrip() for line in output.strip().splitlines())
    return {"path": os.path.realpath(executable), "version": normalized}


@dataclass(frozen=True)
class OutputRequirement:
    key: str
    path: str
    strategy: str = "stat"
    allow_empty: bool = False
    validator: str | None = None
    required_header: tuple[str, ...] = ()

    def normalized_path(self) -> str:
        return _absolute_path(self.path)


def _validate_output(requirement: OutputRequirement) -> str | None:
    path = requirement.normalized_path()
    if not os.path.isfile(path):
        return "output_missing"
    if not requirement.allow_empty and os.path.getsize(path) == 0:
        return "output_empty"

    try:
        if requirement.validator == "json":
            with open(path, "r", encoding="utf-8") as handle:
                json.load(handle)
        elif requirement.validator == "gzip":
            with open(path, "rb") as handle:
                if handle.read(2) != b"\x1f\x8b" or (
                    requirement.strategy == "gzip_crc32"
                    and os.path.getsize(path) < 18
                ):
                    return "output_invalid_gzip"
        elif requirement.validator == "fastq":
            # CRISPResso may emit either gzip-compressed FASTQ or plain FASTQ
            # with a historical .fastq.gz suffix.  Inspect only the format
            # signature so validation remains independent of payload size.
            with open(path, "rb") as handle:
                signature = handle.read(2)
            if signature != b"\x1f\x8b" and not signature.startswith(b"@"):
                return "output_invalid_fastq"
        elif requirement.validator == "bam":
            completed = subprocess.run(
                ["samtools", "quickcheck", path],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=False,
            )
            if completed.returncode != 0:
                return "output_invalid_bam"
        elif requirement.validator == "tsv":
            with open(path, "r", encoding="utf-8") as handle:
                header = tuple(handle.readline().rstrip("\n").split("\t"))
            if requirement.required_header and header[: len(requirement.required_header)] != requirement.required_header:
                return "output_invalid_header"
        elif requirement.validator == "cell_counts":
            with open(path, "r", encoding="utf-8") as handle:
                for line in handle:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) != 2 or not fields[0]:
                        return "output_invalid_cell_counts"
                    int(fields[1])
        elif requirement.validator is not None:
            raise ValueError(f"Unknown output validator: {requirement.validator}")
    except (OSError, UnicodeError, ValueError, json.JSONDecodeError):
        return "output_invalid"
    return None


def _fingerprint_requirement(requirement: OutputRequirement) -> dict[str, object]:
    if requirement.strategy == "sha256":
        fingerprint = small_file_fingerprint(requirement.path)
    elif requirement.strategy == "stat":
        fingerprint = large_file_fingerprint(requirement.path)
    elif requirement.strategy == "gzip_crc32":
        fingerprint = gzip_content_fingerprint(requirement.path)
    else:
        raise ValueError(f"Unknown output fingerprint strategy: {requirement.strategy}")
    return {
        "key": requirement.key,
        **fingerprint,
        "allow_empty": requirement.allow_empty,
        "validator": requirement.validator,
    }


@dataclass
class CacheRecord:
    stage: str
    scope: str
    algorithm_version: int
    dependencies: list[dict[str, object]] = field(default_factory=list)
    inputs: dict[str, object] = field(default_factory=dict)
    parameters: dict[str, object] = field(default_factory=dict)
    tools: dict[str, object] = field(default_factory=dict)
    outputs: list[dict[str, object]] = field(default_factory=list)
    result: dict[str, object] = field(default_factory=lambda: {"status": "completed"})
    producer: dict[str, str] = field(default_factory=dict)
    completed_at: str | None = None
    schema_version: int = CACHE_SCHEMA_VERSION

    @property
    def key_payload(self) -> dict[str, object]:
        return {
            "schema_version": self.schema_version,
            "stage": self.stage,
            "scope": self.scope,
            "algorithm_version": self.algorithm_version,
            "dependencies": self.dependencies,
            "inputs": self.inputs,
            "parameters": self.parameters,
            "tools": self.tools,
        }

    @property
    def cache_key(self) -> str:
        return canonical_digest(self.key_payload)

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_version": self.schema_version,
            "stage": self.stage,
            "scope": self.scope,
            "algorithm_version": self.algorithm_version,
            "cache_key": self.cache_key,
            "producer": self.producer,
            "dependencies": _canonical_value(self.dependencies),
            "inputs": _canonical_value(self.inputs),
            "parameters": _canonical_value(self.parameters),
            "tools": _canonical_value(self.tools),
            "outputs": _canonical_value(self.outputs),
            "result": _canonical_value(self.result),
            "completed_at": self.completed_at,
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "CacheRecord":
        required = {"schema_version", "stage", "scope", "algorithm_version", "cache_key"}
        missing = required - set(payload)
        if missing:
            raise ValueError(f"Cache record is missing fields: {sorted(missing)}")
        record = cls(
            schema_version=int(payload["schema_version"]),
            stage=str(payload["stage"]),
            scope=str(payload["scope"]),
            algorithm_version=int(payload["algorithm_version"]),
            producer=dict(payload.get("producer", {})),
            dependencies=list(payload.get("dependencies", [])),
            inputs=dict(payload.get("inputs", {})),
            parameters=dict(payload.get("parameters", {})),
            tools=dict(payload.get("tools", {})),
            outputs=list(payload.get("outputs", [])),
            result=dict(payload.get("result", {})),
            completed_at=payload.get("completed_at") and str(payload["completed_at"]),
        )
        if payload["cache_key"] != record.cache_key:
            raise ValueError("Cache record key does not match its payload")
        return record


@dataclass(frozen=True)
class CacheDecision:
    status: str
    reasons: tuple[str, ...]
    cache_key: str
    record: CacheRecord | None = None

    @property
    def is_hit(self) -> bool:
        return self.status == "hit"


def atomic_write_json(path: os.PathLike[str] | str, payload: Mapping[str, object]) -> str:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=".tmp", dir=str(destination.parent)
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise
    return str(destination)


def _safe_scope_token(scope: str) -> str:
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", scope).strip("._") or "scope"
    return f"{token}--{hashlib.sha256(scope.encode('utf-8')).hexdigest()}"


class CacheManager:
    def __init__(
        self,
        output_root: str,
        config: CacheConfig | None = None,
        *,
        producer_version: str = "unknown",
        event_callback: Callable[[dict[str, object]], None] | None = None,
    ):
        self.output_root = _absolute_path(output_root)
        self.cache_root = self.output_root + ".cache"
        self.config = config or CacheConfig()
        self.producer_version = producer_version
        self.event_callback = event_callback
        self.events: list[dict[str, object]] = []

    def record_path(self, stage: str, scope: str = "run") -> str:
        if scope == "run":
            return os.path.join(self.cache_root, f"{stage}.json")
        return os.path.join(self.cache_root, stage, _safe_scope_token(scope) + ".json")

    def new_record(
        self,
        stage: str,
        scope: str = "run",
        *,
        algorithm_version: int,
        dependencies: Sequence[Mapping[str, object]] = (),
        inputs: Mapping[str, object] | None = None,
        parameters: Mapping[str, object] | None = None,
        tools: Mapping[str, object] | None = None,
    ) -> CacheRecord:
        return CacheRecord(
            stage=stage,
            scope=scope,
            algorithm_version=algorithm_version,
            dependencies=[dict(item) for item in dependencies],
            inputs=dict(inputs or {}),
            parameters=dict(parameters or {}),
            tools=dict(tools or {}),
            producer={"name": "CRISPRSCope", "version": self.producer_version},
        )

    def _read_record(self, path: str) -> CacheRecord:
        with open(path, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
        record = CacheRecord.from_dict(payload)
        if record.schema_version != CACHE_SCHEMA_VERSION:
            raise ValueError(
                f"Unsupported cache schema {record.schema_version}; expected {CACHE_SCHEMA_VERSION}"
            )
        return record

    def load(self, stage: str, scope: str = "run") -> CacheRecord | None:
        if self.config.mode is CacheMode.DISABLED:
            return None
        try:
            return self._read_record(self.record_path(stage, scope))
        except (FileNotFoundError, OSError, ValueError, TypeError, json.JSONDecodeError):
            return None

    def trusted_stage_records(self, stage: str) -> tuple[CacheRecord, ...]:
        """Return supported per-scope records whose paths match their identities.

        This inventory is intentionally stricter than :meth:`load` because its
        results may be used to authorize removal of obsolete stage-owned files.
        """
        if self.config.mode is CacheMode.DISABLED:
            return ()
        directory = Path(self.cache_root, stage)
        if not directory.is_dir():
            return ()
        records: list[CacheRecord] = []
        for path in sorted(directory.glob("*.json")):
            try:
                record = self._read_record(str(path))
                if record.stage != stage or record.scope == "run":
                    raise ValueError("record stage or scope is inconsistent")
                if Path(self.record_path(record.stage, record.scope)) != path:
                    raise ValueError("record filename does not match its identity")
            except (OSError, ValueError, TypeError, json.JSONDecodeError) as error:
                logging.warning(
                    "Ignoring untrusted cache record during stale cleanup: %s (%s)",
                    path,
                    error,
                )
                continue
            records.append(record)
        return tuple(records)

    def remove_record(self, record: CacheRecord) -> bool:
        """Remove the exact record path derived from a trusted record identity."""
        path = Path(self.record_path(record.stage, record.scope))
        return safe_remove_owned(
            path,
            allowed_root=Path(self.cache_root, record.stage),
            expected_name=path.name,
        )

    def _emit(self, decision: CacheDecision, stage: str, scope: str) -> None:
        event: dict[str, object] = {
            "stage": stage,
            "scope": scope,
            "status": decision.status,
            "reasons": list(decision.reasons),
            "cache_key": decision.cache_key,
        }
        self.events.append(event)
        reason = f" reason={','.join(decision.reasons)}" if decision.reasons else ""
        logging.info("CACHE %s stage=%s scope=%s%s", decision.status.upper(), stage, scope, reason)
        if self.event_callback is not None:
            self.event_callback(dict(event))

    def evaluate(
        self, record: CacheRecord, requirements: Sequence[OutputRequirement]
    ) -> CacheDecision:
        if self.config.mode is CacheMode.DISABLED:
            decision = CacheDecision("disabled", ("cache_disabled",), record.cache_key)
            self._emit(decision, record.stage, record.scope)
            return decision
        if self.config.mode is CacheMode.REFRESH:
            decision = CacheDecision("refresh", ("refresh_requested",), record.cache_key)
            self._emit(decision, record.stage, record.scope)
            return decision

        path = self.record_path(record.stage, record.scope)
        if not os.path.isfile(path):
            decision = CacheDecision("miss", ("record_missing",), record.cache_key)
            self._emit(decision, record.stage, record.scope)
            return decision
        try:
            cached = self._read_record(path)
        except (OSError, ValueError, TypeError, json.JSONDecodeError) as error:
            decision = CacheDecision(
                "invalid", (f"record_unreadable:{type(error).__name__}",), record.cache_key
            )
            self._emit(decision, record.stage, record.scope)
            return decision

        if cached.stage != record.stage or cached.scope != record.scope:
            reasons = ("record_identity_mismatch",)
        elif cached.algorithm_version != record.algorithm_version:
            reasons = ("algorithm_version_changed",)
        elif cached.cache_key != record.cache_key:
            reasons = ("cache_key_changed",)
        else:
            reasons_list: list[str] = []
            cached_outputs = {str(item.get("key")): item for item in cached.outputs}
            expected_keys = {item.key for item in requirements}
            if set(cached_outputs) != expected_keys:
                reasons_list.append("output_inventory_changed")
            for requirement in requirements:
                if reasons_list:
                    break
                cached_output = cached_outputs.get(requirement.key)
                if cached_output is None:
                    reasons_list.append(f"output_missing_record:{requirement.key}")
                    break
                if cached_output.get("path") != requirement.normalized_path():
                    reasons_list.append(f"output_path_changed:{requirement.key}")
                    break
                validation_error = _validate_output(requirement)
                if validation_error:
                    reasons_list.append(f"{validation_error}:{requirement.key}")
                    break
                current = _fingerprint_requirement(requirement)
                if requirement.strategy == "sha256":
                    comparable_keys = ("path", "strategy", "size", "sha256")
                elif requirement.strategy == "gzip_crc32":
                    comparable_keys = (
                        "path", "strategy", "size", "crc32", "uncompressed_size",
                    )
                else:
                    comparable_keys = ("path", "strategy", "size", "mtime_ns")
                if any(cached_output.get(key) != current.get(key) for key in comparable_keys):
                    reasons_list.append(f"output_changed:{requirement.key}")
                    break
            reasons = tuple(reasons_list)

        if reasons:
            decision = CacheDecision("invalid", reasons, record.cache_key, cached)
        else:
            decision = CacheDecision("hit", (), record.cache_key, cached)
        self._emit(decision, record.stage, record.scope)
        return decision

    def commit(
        self,
        record: CacheRecord,
        requirements: Sequence[OutputRequirement],
        *,
        result: Mapping[str, object] | None = None,
    ) -> CacheRecord:
        if self.config.mode is CacheMode.DISABLED:
            return record
        outputs = []
        for requirement in requirements:
            validation_error = _validate_output(requirement)
            if validation_error:
                raise RuntimeError(
                    f"Cannot commit cache record for {record.stage}/{record.scope}: "
                    f"{validation_error} ({requirement.path})"
                )
            outputs.append(_fingerprint_requirement(requirement))
        record.outputs = outputs
        record.result = dict(result or {"status": "completed"})
        record.completed_at = datetime.now(timezone.utc).isoformat()
        atomic_write_json(self.record_path(record.stage, record.scope), record.as_dict())
        return record

    @staticmethod
    def dependency(
        record: CacheRecord,
        output_keys: Iterable[str] | None = None,
        *,
        cache_key: str | None = None,
    ) -> dict[str, object]:
        selected = set(output_keys) if output_keys is not None else None
        outputs = {
            str(item["key"]): {
                key: value
                for key, value in item.items()
                if key in {
                    "path", "strategy", "size", "mtime_ns", "sha256",
                    "crc32", "uncompressed_size",
                }
            }
            for item in record.outputs
            if selected is None or item.get("key") in selected
        }
        return {
            "stage": record.stage,
            "scope": record.scope,
            "cache_key": cache_key or record.cache_key,
            "outputs": outputs,
        }


class OutputRootLock:
    """Non-blocking advisory lock for one output root."""

    def __init__(self, output_root: str):
        self.path = _absolute_path(output_root) + ".cache.lock"
        self._handle = None

    def acquire(self) -> "OutputRootLock":
        Path(self.path).parent.mkdir(parents=True, exist_ok=True)
        handle = open(self.path, "a+", encoding="utf-8")
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            handle.seek(0)
            owner = handle.read().strip() or "owner details unavailable"
            handle.close()
            raise RuntimeError(
                f"Another CRISPRSCope process is using this output root ({owner})"
            ) from error
        except BaseException:
            handle.close()
            raise
        handle.seek(0)
        handle.truncate()
        json.dump(
            {
                "hostname": socket.gethostname(),
                "pid": os.getpid(),
                "started_at": datetime.now(timezone.utc).isoformat(),
            },
            handle,
        )
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
        self._handle = handle
        return self

    def release(self) -> None:
        if self._handle is None:
            return
        try:
            fcntl.flock(self._handle.fileno(), fcntl.LOCK_UN)
        finally:
            self._handle.close()
            self._handle = None

    def __enter__(self) -> "OutputRootLock":
        return self.acquire()

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.release()


def safe_remove_owned(
    path: os.PathLike[str] | str,
    *,
    allowed_root: os.PathLike[str] | str,
    expected_name: str | None = None,
) -> bool:
    """Remove one exact stage-owned path without following it outside its root."""
    candidate = Path(os.path.abspath(path))
    root = Path(os.path.abspath(allowed_root)).resolve()
    resolved = candidate.resolve(strict=False)
    try:
        resolved.relative_to(root)
    except ValueError as error:
        if candidate.is_symlink():
            raise ValueError(
                f"Refusing to remove symlink outside cache-owned root: {candidate}"
            ) from error
        raise ValueError(f"Refusing to remove path outside cache-owned root: {candidate}") from error
    if expected_name is not None and candidate.name != expected_name:
        raise ValueError(f"Refusing to remove unexpected cache-owned name: {candidate.name}")
    if not candidate.exists() and not candidate.is_symlink():
        return False
    if candidate.is_dir() and not candidate.is_symlink():
        shutil.rmtree(candidate)
    else:
        candidate.unlink()
    return True
