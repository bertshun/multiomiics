import os
import sys
import time
import signal
import logging
from pathlib import Path
import pandas as pd
import numpy as np

# Optional imports
try:
    import cptac
except ImportError:
    cptac = None
try:
    import GEOparse
except ImportError:
    GEOparse = None
try:
    from tcia_utils import nbia
except ImportError:
    nbia = None

# --------------------------------------
# ログ設定
# --------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
LOG = logging.getLogger("download_all_safe")

BASE = Path("testdir")
RAW = BASE / "raw"
RAW.mkdir(parents=True, exist_ok=True)

def safe_write_text(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(str(text), encoding="utf-8")

# --------------------------------------
# タイムアウトハンドラ
# --------------------------------------
def timeout_handler(signum, frame):
    raise TimeoutError("Operation timed out")

signal.signal(signal.SIGALRM, timeout_handler)

# --------------------------------------
# GEO Download
# --------------------------------------
def fetch_geo():
    if GEOparse is None:
        LOG.warning("GEOparse not installed.")
        return
    outdir = RAW / "GEO"
    outdir.mkdir(parents=True, exist_ok=True)
    gse_id = "GSE4290"  # Glioblastoma dataset
    try:
        signal.alarm(60)
        LOG.info("Fetching GEO dataset: %s", gse_id)
        gse = GEOparse.get_GEO(geo=gse_id, destdir=str(outdir))
        meta_text = f"GSE title: {gse.get_title()}\nSamples: {len(gse.gsms)}"
        safe_write_text(outdir / f"{gse_id}_meta.txt", meta_text)
        LOG.info("GEO %s saved successfully.", gse_id)
    except TimeoutError:
        LOG.warning("GEO fetch timed out.")
    except Exception as e:
        LOG.exception("GEO fetch error: %s", e)
    finally:
        signal.alarm(0)

# --------------------------------------
# CPTAC Download
# --------------------------------------
def fetch_cptac():
    if cptac is None:
        LOG.warning("cptac package not installed.")
        return
    outdir = RAW / "CPTAC"
    outdir.mkdir(parents=True, exist_ok=True)
    try:
        signal.alarm(60)
        datasets = cptac.list_datasets()
        dataset_list = []
        if isinstance(datasets, pd.DataFrame):
            dataset_list = datasets["Cancer"].dropna().unique().tolist()
        elif isinstance(datasets, (list, tuple)):
            dataset_list = list(datasets)
        else:
            dataset_list = [str(datasets)]

        safe_write_text(outdir / "cptac_available_datasets.txt", "\n".join(dataset_list))
        LOG.info("CPTAC datasets detected: %s", dataset_list)

        # 脳腫瘍（gbm）優先
        chosen = next((ds for ds in dataset_list if "gbm" in ds.lower()), None)
        if not chosen:
            chosen = dataset_list[0] if dataset_list else None
        LOG.info("Chosen CPTAC dataset: %s", chosen)

        if chosen:
            try:
                cptac.download(chosen, redownload=False)
                LOG.info("Downloaded metadata for %s (light mode)", chosen)
            except Exception as e:
                LOG.warning("CPTAC download error for %s: %s", chosen, e)
    except TimeoutError:
        LOG.warning("CPTAC fetch timed out.")
    except Exception as e:
        LOG.exception("CPTAC fetch error: %s", e)
    finally:
        signal.alarm(0)

# --------------------------------------
# TCIA Download
# --------------------------------------
def fetch_tcia():
    if nbia is None:
        LOG.warning("tcia_utils not installed.")
        return
    outdir = RAW / "TCIA"
    outdir.mkdir(parents=True, exist_ok=True)
    try:
        signal.alarm(60)
        LOG.info("Fetching TCIA metadata for collection: TCGA-GBM")
        df = nbia.getSeries(collection="TCGA-GBM")
        df.to_csv(outdir / "tcia_gbm_meta.csv", index=False)
        LOG.info("TCIA metadata saved.")
    except TimeoutError:
        LOG.warning("TCIA fetch timed out.")
    except Exception as e:
        LOG.exception("TCIA fetch error: %s", e)
    finally:
        signal.alarm(0)

# --------------------------------------
# NHANES sample placeholder
# --------------------------------------
def fetch_nhanes():
    outdir = RAW / "NHANES"
    outdir.mkdir(parents=True, exist_ok=True)
    LOG.info("Simulating NHANES blood data fetch (placeholder)")
    df = pd.DataFrame({
        "patient_id": [f"p{i}" for i in range(1, 6)],
        "WBC": np.random.normal(5, 1, 5),
        "RBC": np.random.normal(4.5, 0.5, 5),
        "Hemoglobin": np.random.normal(14, 1, 5)
    })
    df.to_csv(outdir / "nhanes_blood_sample.csv", index=False)
    LOG.info("NHANES sample saved.")

# --------------------------------------
# Run all
# --------------------------------------
def main():
    LOG.info("=== START SAFE DATA FETCH ===")
    fetch_geo()
    fetch_cptac()
    fetch_tcia()
    fetch_nhanes()
    LOG.info("=== COMPLETE ===")

if __name__ == "__main__":
    main()
