import requests
import json
import pandas as pd
import os
import time
from tqdm import tqdm

# =============================
# 設定
# =============================
PROJECT = "TCGA-BRCA"  # 任意のプロジェクト
DATA_TYPE = "Gene Expression Quantification"
SAVE_DIR = "GDC_download"
os.makedirs(SAVE_DIR, exist_ok=True)

# =============================
# 1. メタデータ取得
# =============================
print(f"🔍 {PROJECT} のメタデータを取得中...")

url = "https://api.gdc.cancer.gov/files"
params = {
    "filters": json.dumps({
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [PROJECT]}},
            {"op": "in", "content": {"field": "data_type", "value": [DATA_TYPE]}},
            {"op": "in", "content": {"field": "data_format", "value": ["TXT", "TSV", "CSV", "htseq.counts"]}}
        ]
    }),
    "fields": "file_id,file_name,cases.submitter_id,data_format,data_type,cases.samples.sample_type",
    "format": "JSON",
    "size": "1000"
}

response = requests.get(url, params=params, timeout=60)
response.raise_for_status()
files = response.json()["data"]["hits"]
print(f"✅ {len(files)} 件のファイルが見つかりました。")

# =============================
# 2. 安全ダウンロード関数
# =============================
def safe_download(file_id, save_path, max_retries=5, chunk_size=8192):
    """
    中断時に再試行し、部分ダウンロード再開も可能な安全ダウンローダー
    """
    dl_url = f"https://api.gdc.cancer.gov/data/{file_id}"
    retries = 0

    while retries < max_retries:
        try:
            headers = {}
            downloaded_bytes = os.path.getsize(save_path) if os.path.exists(save_path) else 0

            # ファイルの全体サイズ確認用 HEAD リクエスト
            head = requests.head(dl_url, timeout=30)
            total_size = int(head.headers.get("content-length", 0))

            # ファイルが完全にあるならスキップ
            if downloaded_bytes >= total_size > 0:
                print(f"✅ 既に完全にダウンロード済み: {save_path}")
                return True

            # 再開用ヘッダ
            if downloaded_bytes > 0:
                headers["Range"] = f"bytes={downloaded_bytes}-"

            with requests.get(dl_url, stream=True, headers=headers, timeout=120) as r:
                if r.status_code == 416:
                    print(f"⚠️ {file_id}: ファイルはすでに完全にダウンロード済み (416)。スキップします。")
                    return True

                r.raise_for_status()
                total_to_download = int(r.headers.get("content-length", 0)) + downloaded_bytes

                mode = "ab" if downloaded_bytes > 0 else "wb"
                with open(save_path, mode) as f:
                    pbar = tqdm(
                        total=total_to_download,
                        initial=downloaded_bytes,
                        unit="B",
                        unit_scale=True,
                        desc=os.path.basename(save_path)
                    )

                    for chunk in r.iter_content(chunk_size=chunk_size):
                        if chunk:
                            f.write(chunk)
                            pbar.update(len(chunk))
                    pbar.close()

            return True  # 正常終了

        except requests.exceptions.HTTPError as e:
            if e.response is not None and e.response.status_code == 416:
                print(f"⚠️ {file_id}: 416エラー（ファイル完全）。スキップします。")
                return True
            retries += 1
            print(f"⚠️ HTTPエラー {e}. 再試行 {retries}/{max_retries}")
            time.sleep(5 * retries)

        except (requests.exceptions.ConnectionError,
                requests.exceptions.Timeout,
                requests.exceptions.ChunkedEncodingError) as e:
            retries += 1
            wait = 5 * retries
            print(f"⚠️ 接続中断 ({retries}/{max_retries}) - {e}. {wait}s後に再試行。")
            time.sleep(wait)

    print(f"❌ {file_id} のダウンロードに失敗しました。")
    return False



# =============================
# 3. 一括ダウンロード
# =============================
records = []
for f in files:
    file_id = f["file_id"]
    file_name = f["file_name"]
    sample_id = f["cases"][0]["submitter_id"] if f["cases"] else "Unknown"
    local_path = os.path.join(SAVE_DIR, file_name)

    print(f"\n⬇️ {file_name} をダウンロード中...")
    ok = safe_download(file_id, local_path)

    records.append({
        "file_id": file_id,
        "file_name": file_name,
        "sample_id": sample_id,
        "data_type": f.get("data_type"),
        "data_format": f.get("data_format"),
        "local_path": local_path,
        "status": "success" if ok else "failed"
    })

# =============================
# 4. 結果をDataFrame化・保存
# =============================
df = pd.DataFrame(records)
meta_path = os.path.join(SAVE_DIR, f"{PROJECT}_metadata.csv")
df.to_csv(meta_path, index=False)
print(f"\n💾 メタデータを保存しました → {meta_path}")

# 失敗ファイル一覧表示
failed = df[df["status"] == "failed"]
if not failed.empty:
    print("\n⚠️ ダウンロードに失敗したファイルがあります:")
    print(failed[["file_id", "file_name"]])
