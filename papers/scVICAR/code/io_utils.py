from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    tmp.replace(path)


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_file(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(block_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def checksum_manifest(root: Path, exclude: Iterable[str] = ("SHA256SUMS",)) -> list[tuple[str, str]]:
    blocked = set(exclude)
    rows: list[tuple[str, str]] = []
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        rel = path.relative_to(root).as_posix()
        if rel in blocked:
            continue
        rows.append((sha256_file(path), rel))
    return rows


def write_checksum_manifest(root: Path) -> Path:
    target = root / "SHA256SUMS"
    rows = checksum_manifest(root)
    target.write_text("".join(f"{digest}  {rel}\n" for digest, rel in rows), encoding="utf-8")
    return target


def verify_checksum_manifest(root: Path) -> None:
    manifest = root / "SHA256SUMS"
    if not manifest.exists():
        raise FileNotFoundError(manifest)
    for line in manifest.read_text(encoding="utf-8").splitlines():
        digest, rel = line.split("  ", 1)
        path = root / rel
        if not path.is_file() or sha256_file(path) != digest:
            raise ValueError(f"Checksum mismatch: {path}")


def git_revision(project_root: Path) -> dict[str, Any]:
    def capture(*args: str) -> str:
        result = subprocess.run(
            ["git", "-C", str(project_root), *args],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        return result.stdout.strip()

    return {
        "commit": capture("rev-parse", "HEAD"),
        "branch": capture("branch", "--show-current"),
        "dirty": bool(capture("status", "--porcelain")),
    }


def disk_free_gib(path: Path) -> float:
    return shutil.disk_usage(path).free / (1024**3)


def require_disk_space(path: Path, minimum_gib: float = 5.0) -> None:
    free = disk_free_gib(path)
    if free < minimum_gib:
        raise RuntimeError(f"Only {free:.2f} GiB free at {path}; minimum is {minimum_gib:.2f} GiB")


def safe_environment() -> dict[str, str]:
    keys = ["CUDA_VISIBLE_DEVICES", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"]
    return {key: os.environ[key] for key in keys if key in os.environ}

