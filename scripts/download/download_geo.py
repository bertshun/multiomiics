# scripts/download/download_geo.py
import GEOparse
from pathlib import Path
from scripts.utils import LOG

OUTDIR = Path("data/GEO")
OUTDIR.mkdir(parents=True, exist_ok=True)

def download_gse(gse_id):
    LOG.info(f"Downloading {gse_id}")
    gse = GEOparse.get_GEO(geo=gse_id, destdir=str(OUTDIR))
    # attempt to write expression matrix
    try:
        gse_table = gse.pivot_samples("VALUE")
        gse_table.reset_index(inplace=True)
        gse_table.to_csv(OUTDIR/f"{gse_id}_expr.csv", index=False)
        LOG.info("Saved GEO expression")
    except Exception as e:
        LOG.warning("Could not pivot expression: %s", e)
    return gse

if __name__ == "__main__":
    # user: list GSEs related to GBM/LGG
    gse_list = ["GSE4290","GSE16011"]  # examples; replace as needed
    for g in gse_list:
        download_gse(g)
    LOG.info("GEO downloads done")
