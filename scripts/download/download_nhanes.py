# scripts/download/download_nhanes.py
import requests
from pathlib import Path
from scripts.utils import LOG

OUTDIR = Path("data/EXTERNAL/nhanes")
OUTDIR.mkdir(parents=True, exist_ok=True)

BASE = "https://wwwn.cdc.gov/Nchs/Nhanes/"

def download_cycle(cycle="2017-2018", filecodes=None):
    if filecodes is None:
        filecodes = ["DEMO", "BMX"]  # demographics, body measures
    for f in filecodes:
        fname = f"{f}_{cycle[-1]}.XPT"  # pattern used earlier; validate per file
        url = f"{BASE}{cycle}/{fname}"
        LOG.info(f"Attempt download: {url}")
        r = requests.get(url, timeout=30)
        if r.status_code == 200:
            p = OUTDIR / fname
            with open(p, "wb") as fh:
                fh.write(r.content)
            LOG.info(f"Saved {p}")
        else:
            LOG.warning(f"Not found: {url} status {r.status_code}")

if __name__ == "__main__":
    download_cycle("2017-2018", ["DEMO","BMX","LAB10"])
