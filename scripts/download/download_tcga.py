# scripts/download/download_tcga.py
import requests, json, os
from pathlib import Path
from scripts.utils import LOG, env

GDC_API = "https://api.gdc.cancer.gov"
OUTDIR = Path("data/TCGA")

def query_tcga_cases(projects=["TCGA-GBM","TCGA-LGG"], data_types=None, fields=None, size=10000):
    filters = {
        "op":"in",
        "content":[ {"field":"cases.project.project_id","value": projects} ]
    }
    params = {"filters": json.dumps(filters), "size": size}
    if fields:
        params["fields"] = fields
    r = requests.get(f"{GDC_API}/cases", params=params)
    r.raise_for_status()
    return r.json()

def download_clinical(projects=["TCGA-GBM","TCGA-LGG"]):
    LOG.info("Fetching TCGA clinical via GDC API (cases)")
    # fields to get
    fields = "case_id,submitter_id,diagnoses.age_at_diagnosis,demographic.gender"
    res = query_tcga_cases(projects=projects, fields=fields)
    records = res.get('data', {}).get('hits', [])
    out = OUTDIR / "clinical_gbm_lgg.json"
    with open(out, "w") as f:
        json.dump(records, f)
    LOG.info(f"Written {out}")
    # convert to CSV basic table
    import pandas as pd
    rows = []
    for c in records:
        rid = c.get('submitter_id') or c.get('case_id')
        diag = c.get('diagnoses',[{}])[0]
        rows.append({
            "patient_id": rid,
            "age": diag.get("age_at_diagnosis"),
            "gender": c.get("demographic",{}).get("gender")
        })
    pd.DataFrame(rows).to_csv(OUTDIR/"clinical_gbm_lgg.csv", index=False)
    LOG.info("Saved clinical CSV")

if __name__ == "__main__":
    ensure = Path("data/TCGA")
    ensure.mkdir(parents=True, exist_ok=True)
    download_clinical()
    LOG.info("TCGA metadata fetch complete. For large files use gdc-client with manifest from GDC portal.")
