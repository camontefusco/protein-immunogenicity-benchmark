"""Fast repository integrity checks; does not rerun the publication analyses."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]

def main() -> None:
    required = [
        ROOT / "data/raw/tcell_table_export_1769046013.csv",
        ROOT / "data/curated/assay_level_with_year.csv",
        ROOT / "data/curated/peptide_level_hiconf_with_year.csv",
        ROOT / "data/curated/context_dataset_and_splits.csv",
        ROOT / "data/splits/peptide_split_assignments_v1.csv",
        ROOT / "results/CORRECTED_RESULTS_SUMMARY.md",
    ]
    missing = [str(p.relative_to(ROOT)) for p in required if not p.exists()]
    if missing:
        raise SystemExit("Missing required artifacts: " + ", ".join(missing))
    raw = ROOT / "data/raw/tcell_table_export_1769046013.csv"
    digest = hashlib.sha256(raw.read_bytes()).hexdigest()
    expected = "8f731092a5677180301592500599d1c2fc2c544ffb3426ecd60ca356ebbade02"
    if digest != expected:
        raise SystemExit(f"Raw export checksum mismatch: {digest}")
    config = json.loads((ROOT / "configs/validation_designs.json").read_text())
    if not config:
        raise SystemExit("Validation configuration is empty")
    print("Repository integrity: OK")
    print(f"Raw export SHA-256: {digest}")

if __name__ == "__main__":
    main()
