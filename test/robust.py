import GEOparse
import pandas as pd
from sklearn.preprocessing import RobustScaler

# --- 1. GEO データの取得と整形 ---
# 脳がん関連の GEO シリーズ例を設定
geo_series_list = ["GSE4290", "GSE50161"]  # 必要に応じて追加
geo_dfs = []

for gse_id in geo_series_list:
    print(f"Downloading {gse_id} ...")
    gse = GEOparse.get_GEO(geo=gse_id, destdir="./GEOdata", silent=True)
    
    # サンプルメタデータ取得
    sample_data = gse.phenotype_data.copy()
    
    # 欠損値処理
    sample_data = sample_data.fillna("NA")
    
    # 特定列を整理（例: 病理診断を統一）
    if 'characteristics_ch1.0.Histopathological diagnostic' in sample_data.columns:
        sample_data['histology'] = sample_data['characteristics_ch1.0.Histopathological diagnostic'].str.upper()
    
    geo_dfs.append(sample_data)

# GEO の全シリーズを結合
geo_combined = pd.concat(geo_dfs, axis=0, ignore_index=True)

# --- 2. NHANES データ読み込み ---
nhanes_data = pd.read_csv("./processed_data/NHANES_processed.csv")

# --- 3. Robust Scaling ---
scaler = RobustScaler()

# GEO の数値列のみスケーリング
geo_numeric_cols = geo_combined.select_dtypes(include='number').columns.tolist()
if geo_numeric_cols:
    geo_combined[geo_numeric_cols] = scaler.fit_transform(geo_combined[geo_numeric_cols])

# NHANES の数値列のみスケーリング
nhanes_numeric_cols = nhanes_data.select_dtypes(include='number').columns.tolist()
if nhanes_numeric_cols:
    nhanes_data[nhanes_numeric_cols] = scaler.fit_transform(nhanes_data[nhanes_numeric_cols])

# --- 4. データ保存 ---
geo_combined.to_csv("./processed/GEO_brain_cancer.csv", index=False)
nhanes_data.to_csv("./processed/NHANES_scaled.csv", index=False)

print("GEO と NHANES の整形・スケーリングが完了しました。")
