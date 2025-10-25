# scripts/utils.py
import os
import logging
from pathlib import Path
from dotenv import load_dotenv

load_dotenv()

LOG = logging.getLogger("brain_etl")
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

def ensure_dirs():
    base = Path.cwd()
    for d in ["data/TCGA","data/GEO","data/CPTAC","data/TCIA","data/EXTERNAL","results"]:
        p = base / d
        p.mkdir(parents=True, exist_ok=True)
        LOG.info(f"Ensured {p}")

def env(key, default=None):
    return os.environ.get(key, default)
