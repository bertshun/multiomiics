#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
import numpy as np
import logging

# ログ設定
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
LOG = logging.getLogger("brain_data_merge")

RAW = Path("raw")
PROCESSED = Path("processed")
PROCESSED.mkdir(exist_ok=True)

# 欠損値整形関数
def clean_df(df: pd.DataFrame) -> pd.DataFrame:
    for col in df.columns:
        if pd.api.types.is_numeric_dtype(df[col]):
            df[col] = df[col].fillna(df[col].median())
        else:
            df[col] = df[col].fillna("NA")
    return df

# すべての脳腫瘍データを格納するリスト
all_brain_dfs = []

# -------------------------
# CPTAC 脳腫瘍データ取得
# -------------------------
try:
    import cptac
except ImportError:
    LOG.error("cptac パッケージが見つかりません。pip install cptac でインストールしてください。")
    exit(1)

brain_datasets = ["Gbm", "Lgg"]  # 追加可能

for dataset_name in brain_datasets:
    LOG.info("Fetching CPTAC dataset: %s", dataset_name)
    try:
        dataset_cls = getattr(cptac, dataset_name)
        ds = dataset_cls()
        df_clinical = ds.get_clinical()
        df_clean = clean_df(df_clinical)
        df_clean["cohort"] = f"CPTAC_{dataset_name}"  # cohort カラム追加
        all_brain_dfs.append(df_clean)
    except Exception as e:
        LOG.exception("Failed to fetch or process CPTAC dataset %s: %s", dataset_name, e)

# -------------------------
# NHANES / GEO データ取得
# -------------------------
COHORTS = ["NHANES", "GEO"]

for cohort in COHORTS:
    cohort_dir = RAW / cohort
    csv_files = list(cohort_dir.glob("*.csv"))
    if not csv_files:
        LOG.warning("[%s] No CSV files found in %s", cohort, cohort_dir)
        continue

    LOG.info("[%s] Found %d CSV files", cohort, len(csv_files))
    for csv_file in csv_files:
        try:
            df = pd.read_csv(csv_file)
            df_clean = clean_df(df)
            df_clean["cohort"] = cohort
            all_brain_dfs.append(df_clean)
        except Exception as e:
            LOG.exception("Failed to process %s: %s", csv_file, e)

# -------------------------
# 統合して保存
# -------------------------
if all_brain_dfs:
    merged_df = pd.concat(all_brain_dfs, ignore_index=True, sort=False)
    out_file = PROCESSED / "brain_tumor_combined.csv"
    merged_df.to_csv(out_file, index=False)
    LOG.info("Saved merged brain tumor dataset: %s (rows=%d, cols=%d)", out_file, len(merged_df), merged_df.shape[1])
else:
    LOG.warning("No brain tumor data was found to merge.")
