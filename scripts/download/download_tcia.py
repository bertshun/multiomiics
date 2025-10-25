# scripts/download/download_tcia.py
from tcia_utils import nbia
import pandas as pd
from pathlib import Path
from scripts.utils import LOG, env

OUTDIR = Path("data/TCIA")
OUTDIR.mkdir(parents=True, exist_ok=True)

def fetch_series(collection="TCGA-GBM"):
    LOG.info(f"Fetching TCIA series for {collection}")
    try:
        series = nbia.getSeries(collection=collection)
        df = pd.DataFrame(series)
        df.to_csv(OUTDIR/f"{collection}_series.csv", index=False)
        LOG.info("Saved TCIA series metadata")
    except Exception as e:
        LOG.error("TCIA fetch error: %s", e)

if __name__ == "__main__":
    fetch_series("TCGA-GBM")
