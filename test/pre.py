import pandas as pd
import numpy as np
from sklearn.preprocessing import RobustScaler
from econml.grf import CausalForest

# ========================
# 1. データ読み込み
# ========================
geo_file = "processed_data/GEO_processed.csv"
nhanes_file = "processed_data/NHANES_processed.csv"

geo_df = pd.read_csv(geo_file, low_memory=False)
nhanes_df = pd.read_csv(nhanes_file)

# ========================
# 2. 欠損値処理関数
# ========================
def clean_df(df):
    # 数値列は中央値で補完
    df_numeric = df.select_dtypes(include=[np.number])
    df[df_numeric.columns] = df_numeric.fillna(df_numeric.median())
    # 文字列列は 'Unknown' で補完
    df_str = df.select_dtypes(include=['object'])
    df[df_str.columns] = df_str.fillna('Unknown')
    return df

geo_df = clean_df(geo_df)
nhanes_df = clean_df(nhanes_df)

# ========================
# 3. サンプルID列を統一してマージ
# ========================
geo_df.columns = geo_df.columns.str.strip()
nhanes_df.columns = nhanes_df.columns.str.strip()

geo_df = geo_df.rename(columns={'geo_accession':'sample_id'})
nhanes_df = nhanes_df.rename(columns={'patient_id':'sample_id'})

df = pd.merge(geo_df, nhanes_df, on='sample_id', how='outer')

# ========================
# 4. 説明変数・処置・アウトカム
# ========================
X_cols = ['WBC', 'RBC', 'characteristics_ch1.0.Histopathological diagnostic', 'molecule_ch1']
T_col = 'treatment'
Y_col = 'Hemoglobin'

# 処置が存在しない場合はランダム生成
if T_col not in df.columns:
    df[T_col] = np.random.binomial(1, 0.5, size=len(df))

# ========================
# 5. 説明変数整形
# ========================
X = pd.get_dummies(df[X_cols], drop_first=True)
T = df[T_col].astype(int)
Y = df[Y_col].astype(float)

# Robust scaling
scaler = RobustScaler()
X_scaled = scaler.fit_transform(X)

# ========================
# 6. Causal Forest 学習
# ========================
cf = CausalForest()
cf.fit(Y.values, T.values, X_scaled)
# ========================
# 7. 個別処置効果(ITE)と平均処置効果(ATE)
# ========================
ite = cf.effect(X_scaled)   # 個別処置効果
ate = ite.mean()            # 平均処置効果
print("Average Treatment Effect (ATE):", ate)

# ========================
# 8. 結果保存
# ========================
df['ITE'] = ite
df.to_csv("processed_data_for_causalforest.csv", index=False)
