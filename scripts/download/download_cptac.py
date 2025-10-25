# scripts/download/download_cptac.py
import cptac
from pathlib import Path
from scripts.utils import LOG

OUTDIR = Path("data/CPTAC")
OUTDIR.mkdir(parents=True, exist_ok=True)

def download_brca_or_available():
    LOG.info("Listing CPTAC datasets")
    datasets = cptac.list_datasets()
    LOG.info(f"Available: {datasets}")
    # CPTAC may not have brain; attempt to download brain proteomics if exists, otherwise leave
    try:
        # example: if there is a 'Gbml' dataset – adjust to reality
        # many CPTAC datasets are by tumor types; use whichever is appropriate
        # For demo, we call brca as placeholder
        cptac.download()
        brca = cptac.Brca()
        proteomics = brca.get_proteomics()
        proteomics.to_csv(OUTDIR/"proteomics_brca.csv")
        LOG.info("Saved CPTAC proteomics (BRCA placeholder)")
    except Exception as e:
        LOG.warning("CPTAC download error or dataset not present: %s", e)

if __name__ == "__main__":
    download_brca_or_available()
