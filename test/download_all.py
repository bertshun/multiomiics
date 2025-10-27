#!/usr/bin/env python3
"""
brain_cancer_fetch_and_stats.py

目的:
  - TCGA (public metadata), GEO (GSE), CPTAC (meta), TCIA (series meta),
    PRIDE (projects meta), NHANES (public XPT) から脳腫瘍関連の公開データを取得（APIキー不要）
  - 各コホートごとに簡易クレンジング、RobustScaler 正規化、z-score、binary significance を算出
  - 結果を results/<cohort>_processed.parquet として保存

使い方:
  1) 仮想環境作成
     python -m venv venv
     source venv/bin/activate
     pip install -r requirements.txt
  2) 実行
     python brain_cancer_fetch_and_stats.py

注意:
 - controlled-access データ（dbGaP等）は取得しません。
 - 実データの列名は多様なので、必要に応じて列名マッピングを行ってください。
"""

import os
import sys
import time
import json
import logging
from pathlib import Path
from typing import List, Tuple, Dict

import requests
import pandas as pd
import numpy as np
from sklearn.preprocessing import RobustScaler
from scipy import stats
import time, urllib.parse

# Optional 3rd-party libs: GEOparse, cptac, tcia_utils, pyreadstat
try:
    import GEOparse
except Exception:
    GEOparse = None
try:
    import cptac
except Exception:
    cptac = None
try:
    from tcia_utils import nbia
except Exception:
    nbia = None
try:
    import pyreadstat
except Exception:
    pyreadstat = None

# --- config paths
ROOT = Path.cwd()
RAW = ROOT / "testdir" / "raw"
RESULTS = ROOT / "results"
for p in (RAW, RESULTS):
    p.mkdir(parents=True, exist_ok=True)

# Logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
LOG = logging.getLogger("brain_etl")

# -------------------------
# Helper utilities
# -------------------------
def safe_write_text(path: Path, s: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(s, encoding="utf-8")
    LOG.info("Wrote %s", path)

def safe_write_bytes(path: Path, b: bytes):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as fh:
        fh.write(b)
    LOG.info("Wrote %s (%d bytes)", path, path.stat().st_size)

######### 修正版：TCGA の expression ファイルを探索して小さな TSV/CSV を落とす #########

def fetch_tcga_expression(projects=["TCGA-GBM","TCGA-LGG"], outdir=RAW/"TCGA_expr"):
    outdir.mkdir(parents=True, exist_ok=True)
    LOG.info("Searching TCGA files for expression data: %s", projects)
    base_files = "https://api.gdc.cancer.gov/files"
    # Query for files of type 'Gene Expression Quantification' (open-access)
    filters = {
        "op":"in",
        "content":[{"field":"cases.project.project_id","value":projects},
                   {"field":"files.data_type","value":["Gene Expression Quantification","Transcriptome Profiling"]}]
    }
    params = {"filters": json.dumps(filters), "size": 1000, "fields":"file_id,file_name,access,md5sum,cases.samples.submitter_id,cases.submitter_id"}
    try:
        r = requests.get(base_files, params=params, timeout=60)
        r.raise_for_status()
        j = r.json()
        hits = j.get("data", {}).get("hits", [])
        LOG.info("TCGA files found: %d", len(hits))
        # For each hit, try to find a downloadable URL (data repository or gdc API)
        for h in hits:
            fname = h.get("file_name")
            file_id = h.get("file_id")
            access = h.get("access","")
            # Skip controlled access
            if access.lower()=="controlled":
                LOG.info("Skipping controlled file %s", fname)
                continue
            # Try to download via GDC data endpoint (small files only)
            dl_url = f"https://api.gdc.cancer.gov/data/{file_id}"
            try:
                r2 = requests.get(dl_url, timeout=60, stream=True)
                if r2.status_code == 200 and int(r2.headers.get("Content-Length", "0")) < 100*1024*1024:
                    # save streaming to file (protect memory)
                    outp = outdir / fname
                    with open(outp, "wb") as fh:
                        for chunk in r2.iter_content(1024*1024):
                            fh.write(chunk)
                    LOG.info("Downloaded TCGA file %s", outp)
                else:
                    LOG.info("TCGA file %s too large or unavailable (status=%s, len=%s)", fname, r2.status_code, r2.headers.get("Content-Length"))
            except Exception as e:
                LOG.warning("Failed to download TCGA file %s: %s", fname, e)
        # List downloaded csv/tsv
        downloaded = list(outdir.glob("*"))
        LOG.info("TCGA expr downloaded files: %d", len(downloaded))
    except Exception as e:
        LOG.exception("TCGA files query failed: %s", e)
    return outdir


# -------------------------
# 3) CPTAC (metadata only, no heavy download)
# -------------------------
def fetch_cptac_meta() -> Path:
    outdir = RAW / "CPTAC"
    outdir.mkdir(parents=True, exist_ok=True)
    if cptac is None:
        LOG.warning("cptac package not installed. Skipping CPTAC fetch.")
        return outdir
    try:
        datasets = cptac.list_datasets()
        # datasets may be list or DataFrame
        if hasattr(datasets, "to_dict"):
            # DataFrame-like
            try:
                ds = list(datasets["Cancer"].dropna().unique())
            except Exception:
                ds = list(datasets.columns)
        elif isinstance(datasets, (list, tuple)):
            ds = list(datasets)
        else:
            ds = [str(datasets)]
        safe_write_text(outdir / "cptac_datasets.txt", "\n".join(map(str, ds)))
        LOG.info("CPTAC datasets: %s", ds)
        # avoid heavy ds_obj instantiation; optionally call metadata
    except Exception as e:
        LOG.exception("CPTAC fetch error: %s", e)
    return outdir

# -------------------------
# 4) TCIA series metadata fetch (public)
# -------------------------
# --- 追加修正版 fetch_tcia_series（タイムアウト短縮＋安全化） ---
######### 修正版：TCIA - 取得はメタのみでOK。TCIA からは DICOM を落として radiomics を自前で抽出する方針 #########
def fetch_tcia_safe(collection="TCGA-GBM", outdir=RAW/"TCIA"):
    outdir.mkdir(parents=True, exist_ok=True)
    base = "https://services.cancerimagingarchive.net/services/v4/TCIA/query/getSeries"
    try:
        r = requests.get(base, params={"Collection": collection}, timeout=10)
        if r.status_code == 200:
            text = r.text
            # Save truncated xml/json
            if len(text) > 200000:
                text = text[:200000] + "\n<!-- truncated -->"
            safe_write_text(outdir / f"{collection}_series_trunc.xml", text)
            LOG.info("Saved TCIA truncated xml")
        else:
            LOG.warning("TCIA status %s", r.status_code)
    except requests.exceptions.Timeout:
        LOG.warning("TCIA timed out")
    except Exception as e:
        LOG.exception("TCIA error: %s", e)
    # If you need DICOMs: use TCIA NBIA NBIAClient or HTTPS FTP, requires more code & large storage
    return outdir



# --- スクリプト最後に追記（結果確認用） ---
def show_csv_summary():
    """Fetch完了後、各 cohort の .csv を pandas で開いて df.head() を表示"""
    import pandas as pd

    cohort_dirs = {
        "TCGA": RAW / "TCGA",
        "GEO": RAW / "GEO",
        "CPTAC": RAW / "CPTAC",
        "TCIA": RAW / "TCIA",
        "PRIDE": RAW / "PRIDE",
        "NHANES": RAW / "NHANES",
    }

    for name, d in cohort_dirs.items():
        csvs = list(d.glob("*.csv"))
        if not csvs:
            LOG.info("[%s] No CSV files found", name)
            continue
        LOG.info("[%s] Found %d CSV files", name, len(csvs))
        for csv in csvs:
            try:
                df = pd.read_csv(csv)
                print(f"\n--- {name}: {csv.name} ---")
                print(df.head(5))
                print(f"(rows={len(df)}, cols={df.shape[1]})")
            except Exception as e:
                LOG.warning("Failed to read %s: %s", csv, e)

#

# -------------------------
# 5) PRIDE projects meta for glioma searches
# -------------------------
######### 修正版：PRIDE - プロジェクトメタから「実データのftpリンク」を抽出して小さなファイルを落とす #########
def fetch_pride_and_files(query="glioblastoma", outdir=RAW/"PRIDE_files", max_projects=10):
    outdir.mkdir(parents=True, exist_ok=True)
    LOG.info("PRIDE search for: %s", query)
    api = "https://www.ebi.ac.uk:443/pride/ws/archive/v2/projects"
    params = {"pageSize": 50, "page": 0, "q": query}
    projects = []
    try:
        while True:
            r = requests.get(api, params=params, timeout=30)
            if r.status_code != 200:
                LOG.warning("PRIDE API returned %s", r.status_code)
                break
            data = r.json()
            items = data.get("_embedded", {}).get("projects", [])
            projects.extend(items)
            if data.get("_links", {}).get("next") is None or params["page"] > 10:
                break
            params["page"] += 1
            time.sleep(0.1)
    except Exception as e:
        LOG.exception("PRIDE search failed: %s", e)
    LOG.info("PRIDE projects found: %d", len(projects))
    # iterate top-N projects, get files for each
    count = 0
    for p in projects:
        if count >= max_projects:
            break
        acc = p.get("accession") or p.get("projectAccession")
        if not acc:
            continue
        count += 1
        try:
            files_api = f"https://www.ebi.ac.uk/pride/ws/archive/file/listProjectFiles/{acc}"
            r2 = requests.get(files_api, timeout=30)
            if r2.status_code != 200:
                LOG.warning("PRIDE files API %s returned %s", files_api, r2.status_code)
                continue
            file_items = r2.json().get("list", [])
            # For each file, if ftpDownloadLink exists and it's small (txt/csv/tsv), download
            for fi in file_items:
                ftp = fi.get("ftpDownloadLink") or fi.get("downloadLink")
                fname = fi.get("fileName") or urllib.parse.urlsplit(ftp).path.split("/")[-1] if ftp else None
                if not ftp or not fname:
                    continue
                # prefer small text-like files
                if fname.lower().endswith((".txt",".csv",".tsv",".mzid",".mzML",".gz")):
                    outp = outdir / f"{acc}__{fname}"
                    try:
                        r3 = requests.get(ftp, timeout=60, stream=True)
                        if r3.status_code == 200:
                            # if too large (>200MB) skip
                            clen = int(r3.headers.get("Content-Length","0"))
                            if clen > 200*1024*1024:
                                LOG.info("Skipping large file %s (size=%s)", fname, clen)
                                continue
                            with open(outp, "wb") as fh:
                                for chunk in r3.iter_content(1024*1024):
                                    fh.write(chunk)
                            LOG.info("Downloaded PRIDE file %s", outp)
                        else:
                            LOG.warning("PRIDE file %s unreachable (status=%s)", ftp, r3.status_code)
                    except Exception as e:
                        LOG.warning("Failed download %s: %s", ftp, e)
        except Exception as e:
            LOG.exception("Failed fetching files for project %s: %s", acc, e)
    return outdir

# -------------------------
# 7) Simple ETL/statistics processing function
#    - Input: directory with CSV files or table-like files
#    - Output: processed parquet with scaled numeric cols, zscore, binary_signif
# -------------------------
def process_cohort_dir(cohort_name: str, cohort_dir: Path, out_dir: Path = RESULTS,
                       numeric_only: bool = False, chunk_rows: int = 100000):
    """
    Process CSV-like tables found in cohort_dir.
    For simplicity: find all CSV files in cohort_dir, concat (careful with memory),
    select numeric columns, apply RobustScaler, compute zscore (per column),
    create binary_signif column suffix _sig (|z|>2).
    Save results to results/<cohort_name>_processed.parquet
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_files = list(cohort_dir.glob("*.csv"))
    if not csv_files:
        LOG.warning("No CSV files found for cohort %s in %s", cohort_name, cohort_dir)
        return None
    LOG.info("Processing cohort %s with %d CSV files", cohort_name, len(csv_files))
    # read in chunks to avoid memory explosion: accumulate numeric summary then transform
    # Strategy: compute column-wise median and IQR across files by streaming, then scale per file
    numeric_cols = set()
    col_medians = {}
    col_iqr = {}

    # First pass: discover numeric columns and collect sample statistics using downsampling
    sample_frames = []
    for f in csv_files:
        try:
            # read small sample rows to infer dtypes
            df_sample = pd.read_csv(f, nrows=1000)
            sample_frames.append(df_sample)
            for c in df_sample.columns:
                if pd.api.types.is_numeric_dtype(df_sample[c]):
                    numeric_cols.add(c)
        except Exception as e:
            LOG.warning("Skipping sample read for %s: %s", f, e)

    numeric_cols = sorted(list(numeric_cols))
    if not numeric_cols:
        LOG.warning("No numeric columns found for %s. Skipping.", cohort_name)
        return None
    LOG.info("Identified numeric columns: %s", numeric_cols[:10])

    # Compute robust stats (median, IQR) across files by sampling blocks
    medians_accum = {c: [] for c in numeric_cols}
    q1_accum = {c: [] for c in numeric_cols}
    q3_accum = {c: [] for c in numeric_cols}

    for f in csv_files:
        LOG.info("Sampling for stats from %s", f)
        try:
            for chunk in pd.read_csv(f, usecols=numeric_cols, chunksize=50000):
                # downsample chunk to limit memory
                if len(chunk) > 5000:
                    chunk = chunk.sample(2000, random_state=0)
                for c in numeric_cols:
                    col = chunk[c].dropna()
                    if not col.empty:
                        medians_accum[c].append(col.median())
                        q1_accum[c].append(col.quantile(0.25))
                        q3_accum[c].append(col.quantile(0.75))
                break  # only one chunk per file for initial stats
        except Exception as e:
            LOG.warning("Stat sampling failed for %s: %s", f, e)

    for c in numeric_cols:
        med = np.median(medians_accum[c]) if medians_accum[c] else 0.0
        q1 = np.median(q1_accum[c]) if q1_accum[c] else 0.0
        q3 = np.median(q3_accum[c]) if q3_accum[c] else 0.0
        iqr = q3 - q1 if (q3 - q1) != 0 else 1.0
        col_medians[c] = med
        col_iqr[c] = iqr

    LOG.info("Computed medians and IQRs for %d cols", len(col_medians))

    # Second pass: transform per-file and append to an output parquet (using pyarrow engine)
    out_path = out_dir / f"{cohort_name}_processed.parquet"
    # We'll accumulate into a list of smaller dataframes then concat write (careful with memory).
    processed_chunks = []
    for f in csv_files:
        LOG.info("Transforming file %s", f)
        try:
            for chunk in pd.read_csv(f, chunksize=chunk_rows):
                # select relevant columns
                numeric_chunk = chunk[numeric_cols].copy()
                # fill missing with median
                for c in numeric_cols:
                    numeric_chunk[c] = numeric_chunk[c].fillna(col_medians.get(c, 0.0))
                # Robust scaling: (x - median) / IQR
                for c in numeric_cols:
                    numeric_chunk[c] = (numeric_chunk[c] - col_medians[c]) / col_iqr[c]
                # compute zscore (standardization)
                z_chunk = numeric_chunk.apply(lambda col: stats.zscore(col, nan_policy='omit'))
                # add suffix columns
                for c in numeric_cols:
                    zc = z_chunk[c]
                    chunk[f"{c}_z"] = zc
                    chunk[f"{c}_sig"] = zc.abs() > 2.0
                # append processed chunk (drop original large non-numeric columns to save memory)
                # Keep identifier columns if present
                id_cols = [c for c in chunk.columns if "id" in c.lower() or c.lower() in ("patient_id","sample_id")]
                keep_cols = id_cols + [col for c in numeric_cols for col in (c, f"{c}_z", f"{c}_sig")]
                keep_cols = [c for c in keep_cols if c in chunk.columns]
                processed_chunks.append(chunk[keep_cols].reset_index(drop=True))
                # flush to disk temp if processed_chunks grows big
                if sum(len(df) for df in processed_chunks) > 500000:
                    # write intermediate and clear
                    tmp_path = out_dir / f"{cohort_name}_part_{int(time.time())}.parquet"
                    pd.concat(processed_chunks, ignore_index=True).to_parquet(tmp_path, index=False)
                    LOG.info("Flushed intermediate parquet %s", tmp_path)
                    processed_chunks = []
        except Exception as e:
            LOG.exception("Processing chunk failed for %s: %s", f, e)

    # Final merge of all processed chunks
    if processed_chunks:
        final_df = pd.concat(processed_chunks, ignore_index=True)
        final_df.to_parquet(out_path, index=False)
        LOG.info("Saved processed cohort parquet: %s (rows=%d cols=%d)", out_path, len(final_df), final_df.shape[1])
    else:
        LOG.warning("No processed data chunks produced for %s", cohort_name)
        final_df = None
    return out_path

# -------------------------
# Main runner orchestrating fetch + process
# -------------------------
def main():
    LOG.info("=== START: Fetch public brain-cancer datasets (no keys) and process ===")
    # 1) TCGA metadata
    try:
        fetch_tcga_metadata()
    except Exception as e:
        LOG.exception("TCGA fetch failed: %s", e)

    # 2) GEO
    try:
        fetch_geo_gse("GSE4290")  # example GBM dataset — adjust as needed
    except Exception as e:
        LOG.exception("GEO fetch failed: %s", e)

    # 3) CPTAC
    try:
        fetch_cptac_meta()
    except Exception as e:
        LOG.exception("CPTAC fetch failed: %s", e)

    # 4) TCIA
    try:
        fetch_tcia_series("TCGA-GBM")
    except Exception as e:
        LOG.exception("TCIA fetch failed: %s", e)

    # 5) PRIDE
    try:
        fetch_pride_projects("glioblastoma OR glioma")
    except Exception as e:
        LOG.exception("PRIDE fetch failed: %s", e)

    # 6) NHANES sample
    try:
        fetch_nhanes_sample("2017-2018", ["DEMO","BMX"])
    except Exception as e:
        LOG.exception("NHANES fetch failed: %s", e)

    # 7) Process each cohort where CSVs exist
    # Example mapping: TCGA -> testdir/raw/TCGA (we created a simple CSV earlier)
    cohort_dirs = {
        "tcga": RAW / "TCGA",
        "geo_gse4290": RAW / "GEO",
        "cptac": RAW / "CPTAC",
        "tcia": RAW / "TCIA",
        "pride": RAW / "PRIDE",
        "nhanes": RAW / "NHANES"
    }

    for cohort, dpath in cohort_dirs.items():
        try:
            process_cohort_dir(cohort, dpath)
        except Exception as e:
            LOG.exception("Processing failed for cohort %s: %s", cohort, e)

#  --- main の最後に追加 ---
    LOG.info("=== COMPLETE ===")
    show_csv_summary()

if __name__ == "__main__":
    main()
