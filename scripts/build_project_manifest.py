from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "09_provenance" / "manifests" / "canonical_project_manifest.csv"


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


files = sorted(
    path for path in ROOT.rglob("*")
    if path.is_file() and path != OUTPUT and "__pycache__" not in path.parts
)
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
with OUTPUT.open("w", newline="", encoding="utf-8") as stream:
    writer = csv.DictWriter(stream, fieldnames=["relative_path", "bytes", "sha256"])
    writer.writeheader()
    for path in files:
        writer.writerow(
            {
                "relative_path": str(path.relative_to(ROOT)),
                "bytes": path.stat().st_size,
                "sha256": digest(path),
            }
        )
print(f"Registered {len(files)} files in {OUTPUT}")
