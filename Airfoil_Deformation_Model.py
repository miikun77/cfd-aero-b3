import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft

# ファイルの読み込み
file_path = r"C:\CFD2\airfoil.dat"  # 適宜ファイルパスを指定してください

# データを読み込んでDataFrameに格納
df = pd.read_csv(file_path, delim_whitespace=True, header=None, names=["x", "y"])

# NaNを含む行を削除
df_clean = df.dropna()

# x, yデータとして抽出
x = df_clean["x"].values
y = df_clean["y"].values

# x=0の位置を検出
zero_index = np.where(x == 0)[0][0]

# x=0を基準に分割
df_lower = df_clean.iloc[:zero_index + 1]  # x=0より前のデータ（下面）
df_upper = df_clean.iloc[zero_index:]       # x=0より後のデータ（上面）

# x, y のデータとして分離
x_upper, y_upper = df_upper["x"].values, df_upper["y"].values
x_lower, y_lower = df_lower["x"].values, df_lower["y"].values

# フーリエ変換
y_upper_fft = fft(y_upper)
y_lower_fft = fft(y_lower)

# #翼圧変形
# y_upper_fft *= 4/5
# y_lower_fft *= 4/5

# 逆フーリエ変換で再構成
y_upper_reconstructed = ifft(y_upper_fft).real
y_lower_reconstructed = ifft(y_lower_fft).real

# 最大翼圧取得
t = y_upper_reconstructed - y_lower_reconstructed
t_max = max(t)

# キャンバー変更（たわみ再現？）
c = y_upper_reconstructed + y_lower_reconstructed
c_max = max(c)


# たわみ再現（テキトー）
y_upper_reconstructed -= x_upper*t_max/10 
y_lower_reconstructed -= x_lower*t_max/10 

# プロット
plt.figure(figsize=(10, 5))
plt.plot(x_upper, y_upper, 'o-', label='Upper Surface Original', linewidth=1)
plt.plot(x_upper, y_upper_reconstructed, 'o--', label='Upper Surface Reconstructed', linewidth=1)
plt.plot(x_lower, y_lower, 'o-', label='Lower Surface Original', linewidth=1)
plt.plot(x_lower, y_lower_reconstructed, '--', label='Lower Surface Reconstructed', linewidth=1)
plt.title('Fourier Series Reconstruction of Airfoil Upper and Lower Surfaces')
plt.xlabel('Chordwise Position (x)')
plt.ylabel('Vertical Position (y)')
plt.legend()
plt.grid(True)
plt.axis('equal')
plt.show()
