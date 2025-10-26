import os
import requests

DATA_SOURCES = {
    "TCGA_GBM": "https://tcga-data.nci.nih.gov/tcga_gbm_clinical.csv",
    "TCGA_LGG": "https://tcga-data.nci.nih.gov/tcga_lgg_clinical.csv",
    "GEO": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE/GSEXXXX/matrix.csv.gz",
    "CPTAC": "https://cptac-data-portal.georgetown.edu/brain_proteomics.csv",
    "TCIA": "https://tcia.at/download/tcia_gbm_lgg_features.csv"
}

os.makedirs("testdir/raw", exist_ok=True)

def download_file(url, dest):
    print(f"Downloading {url}")
    r = requests.get(url)
    if r.status_code == 200:
        with open(dest, "wb") as f:
            f.write(r.content)
    else:
        print(f"Failed: {url}, status={r.status_code}")

for name, url in DATA_SOURCES.items():
    dest_path = f"testdir/raw/{name}.csv"
    download_file(url, dest_path)
    print(f"✅ Saved {dest_path}")
