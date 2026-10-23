from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

DATA_DIR = PROJECT_ROOT / "data"

RAW_DIR = DATA_DIR / "raw"
RESULTS_DIR = DATA_DIR / "results"
REFERENCE_DIR = DATA_DIR / "reference"
EXTERNAL_DIR = DATA_DIR / "external"
MANUAL_DIR = DATA_DIR / "manual"
WORK_DIR = DATA_DIR / "work"
PROCESSED_DIR = DATA_DIR / "processed"

MANIFESTS_DIR = PROJECT_ROOT / "manifests"
LLM_EXPORT_DIR = PROJECT_ROOT / "llm_export"
NOTEBOOKS_DIR = PROJECT_ROOT / "notebooks"
