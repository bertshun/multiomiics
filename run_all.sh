#!/bin/bash
set -e
export PYTHONUNBUFFERED=1

# Ensure .env loaded by utils (dotenv)
python - <<'PY'
from scripts.utils import ensure_dirs
ensure_dirs()
print("Dirs ensured")
PY

echo "Downloading TCGA metadata..."
python scripts/download/download_tcga.py

echo "Downloading GEO datasets..."
python scripts/download/download_geo.py

echo "Downloading CPTAC (if available)..."
python scripts/download/download_cptac.py

echo "Downloading TCIA metadata..."
python scripts/download/download_tcia.py

echo "Downloading PRIDE metadata..."
python scripts/download/download_pride.py

echo "Downloading NHANES cycles..."
python scripts/download/download_nhanes.py

echo "Running ETL..."
python scripts/etl/etl_brain.py

echo "Done. Results in results/brain_cancer_etl.csv"
