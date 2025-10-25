# scripts/etl/etl_brain.py
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.preprocessing import RobustScaler
from scripts.utils import LOG, ensure_dirs

ensure_dirs()
OUT = Path("results")
OUT.mkdir(parents=True, exist_ok=True)

def read_if_exists(path):
    p = Path(path)
    if p.exists():
        try:
            if p.suffix.lower() in [".csv",".txt"]:
                return pd.read_csv(p)
            elif p.suffix.lower() in [".json"]:
                return pd.read_json(p)
            else:
                # try reading with pandas generic
                return pd.read_csv(p)
        except Exception as e:
            LOG.warning("Read failed for %s: %s", p, e)
            return pd.DataFrame()
    else:
        LOG.info("File not found: %s", p)
        return pd.DataFrame()

# 1. TCGA
tcga_clinical = read_if_exists("data/TCGA/clinical_gbm_lgg.csv")
tcga_genomic  = read_if_exists("data/TCGA/genomic_gbm_lgg.csv")
if not tcga_genomic.empty:
    tcga_df = tcga_clinical.merge(tcga_genomic, on="patient_id", how="inner")
else:
    tcga_df = tcga_clinical

# 2. GEO
geo_expr = read_if_exists("data/GEO/GSE_expr.csv")  # adjust path/name
if not geo_expr.empty:
    merged = tcga_df.merge(geo_expr, on="patient_id", how="left")
else:
    merged = tcga_df

# 3. CPTAC / PRIDE proteomics
cptac_df = read_if_exists("data/CPTAC/brain_proteomics.csv")
pride_meta = read_if_exists("data/EXTERNAL/pride/pride_projects_meta.csv")
# if proteomics measurements exist, join a few proteins (names may vary)
proteins = ["P53","TP53","VEGFA","IL6"]
if not cptac_df.empty:
    # attempt to map protein columns
    prot_cols = [c for c in cptac_df.columns if any(p.lower() in c.lower() for p in proteins)]
    prot_cols = prot_cols[:10]
    prot_select = ["patient_id"] + prot_cols if prot_cols else ["patient_id"]
    cptac_small = cptac_df.loc[:, [c for c in prot_select if c in cptac_df.columns]]
    cptac_small = cptac_small.rename(columns={c:c.replace("TP53","prot_P53").replace("P53","prot_P53") for c in prot_select})
    merged = merged.merge(cptac_small, on="patient_id", how="left")

# 4. TCIA features
tcia_df = read_if_exists("data/TCIA/TCGA-GBM_series.csv")
if not tcia_df.empty:
    # make sure feature names are consistent
    # assume columns tumor_ratio, necrosis_ratio, inflammation exist
    merged = merged.merge(tcia_df, on="patient_id", how="left")

# 5. NHANES / blood
nhanes = read_if_exists("data/EXTERNAL/nhanes/DEMO_J.XPT")  # if XPT converted to csv earlier
# if using XPT, convert via pyreadstat; here assume CSV exists
if not nhanes.empty:
    merged = merged.merge(nhanes, on="patient_id", how="left")

# 6. Numeric normalization
numeric_cols = [c for c in merged.columns if merged[c].dtype.kind in 'fi']
# Filter to important numeric features if too many
keep = ["age","prot_P53","VEGFA","IL6","tumor_ratio","necrosis_ratio","inflammation","WBC","RBC","hemoglobin"]
numeric_cols = [c for c in keep if c in merged.columns]

if numeric_cols:
    LOG.info("Numeric cols for scaling: %s", numeric_cols)
    scaler = RobustScaler()
    merged[numeric_cols] = merged[numeric_cols].fillna(0)  # temporary fill before scaling
    merged[numeric_cols] = scaler.fit_transform(merged[numeric_cols])
    # now reintroduce NaN where appropriate (but we filled above). Optionally do median-impute:
    for c in numeric_cols:
        if merged[c].isna().sum()>0:
            merged[c].fillna(merged[c].median(), inplace=True)

# 7. Final save
merged.to_csv(OUT/"brain_cancer_etl.csv", index=False)
LOG.info("ETL saved to %s", OUT/"brain_cancer_etl.csv")
