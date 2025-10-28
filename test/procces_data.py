import os
import pandas as pd
from sklearn.preprocessing import RobustScaler

# 出力先ディレクトリ
output_dir = "processed_data"
os.makedirs(output_dir, exist_ok=True)

# データファイル一覧（例: すでにダウンロード済みの CSV を想定）
data_files = {
    "GEO": "processed/GEO/GSE4290_pheno.csv",
    "NHANES": "processed/NHANES/nhanes_blood_sample.csv"
}

# 保存用辞書
cleaned_data = {}

for name, file_path in data_files.items():
    print(f"Processing {name} ...")
    
    # CSV 読み込み
    df = pd.read_csv(file_path)
    
    # 欠損値処理: NaN を "NA" に置換
    df = df.fillna("NA")
    
    # もし数値カラムだけ Robust Scaling したければ別途
    num_cols = df.select_dtypes(include=["int64", "float64"]).columns
    if len(num_cols) > 0:
        scaler = RobustScaler()
        df[num_cols] = scaler.fit_transform(df[num_cols])
    
    # 整形後 CSV 保存
    output_file = os.path.join(output_dir, f"{name}_processed.csv")
    df.to_csv(output_file, index=False)
    print(f"{name} saved to {output_file}")
    
    cleaned_data[name] = df

print("All GEO and NHANES data processed and scaled.")
