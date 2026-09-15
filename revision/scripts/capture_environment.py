from __future__ import annotations

import json
import platform
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "09_provenance" / "environments"
OUT.mkdir(parents=True, exist_ok=True)

freeze = subprocess.run(
    [sys.executable, "-m", "pip", "freeze"],
    check=True,
    capture_output=True,
    text=True,
).stdout
(OUT / "requirements_observed.txt").write_text(freeze, encoding="utf-8")
(OUT / "runtime.json").write_text(
    json.dumps(
        {
            "python_executable": sys.executable,
            "python_version": sys.version,
            "platform": platform.platform(),
            "machine": platform.machine(),
        },
        indent=2,
        sort_keys=True,
    ),
    encoding="utf-8",
)
print(f"Captured environment from {sys.executable}")
