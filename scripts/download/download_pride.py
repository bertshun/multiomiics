# scripts/download/download_pride.py
import requests, pandas as pd
from pathlib import Path
from scripts.utils import LOG, env

OUTDIR = Path("data/EXTERNAL/pride")
OUTDIR.mkdir(parents=True, exist_ok=True)

PRIDE_API = "https://www.ebi.ac.uk/pride/ws/archive/v2/projects"

def fetch_pride_projects(query="cancer"):
    LOG.info("Fetching PRIDE projects list (filtered)")
    params = {"pageSize":100, "page":0, "speciesFilter":"Homo sapiens", "q": query}
    projects = []
    while True:
        r = requests.get(PRIDE_API, params=params, timeout=30)
        r.raise_for_status()
        data = r.json()
        items = data.get("_embedded", {}).get("projects", [])
        projects.extend(items)
        if data.get("_links", {}).get("next") is None:
            break
        params["page"] += 1
    df = pd.DataFrame(projects)
    df.to_csv(OUTDIR/"pride_projects_meta.csv", index=False)
    LOG.info("Saved pride meta")
    return df

if __name__ == "__main__":
    fetch_pride_projects("glioblastoma OR glioma")
