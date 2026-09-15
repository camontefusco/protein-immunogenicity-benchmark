from pathlib import Path
import subprocess
import sys

def test_repository_integrity_script():
    root = Path(__file__).resolve().parents[1]
    result = subprocess.run([sys.executable, str(root / "scripts/audit_repository.py")], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr + result.stdout

def test_streamlit_app_compiles():
    root = Path(__file__).resolve().parents[1]
    result = subprocess.run([sys.executable, "-m", "py_compile", str(root / "apps/iedb_streamlit/app.py")], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
