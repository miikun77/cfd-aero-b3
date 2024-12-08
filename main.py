import subprocess
import os
import re
import csv
import shutil
import numpy as np

def create_directories():
    os.makedirs("build", exist_ok=True)
    os.makedirs("work", exist_ok=True)

def add_snow_effect(parameter):
    # airfoil.datを読み込み
    with open("C:/CFD2/airfoil.dat", "r") as file:
        lines = file.readlines()
    
    # ヘッダーを保持
    header = lines[:2]
    data_lines = lines[2:]
    
    # 座標データを取得
    data = []
    for line in data_lines:
        x_str, y_str = line.strip().split()
        x = float(x_str)
        y = float(y_str)
        data.append([x, y])
    
    # 三次関数を定義 (f(x) = a * x * (1 - x)^2)
    def cubic_function(x, a):
        return a * x * (1 - x)**2
    
    # 上側の翼に雪を追加
    half = len(data) // 2
    for i in range(half, len(data)):
        x, y = data[i]
        data[i][1] += cubic_function(x, parameter)
    
    # データを書き込み
    with open("C:/CFD2/work/airfoil-snowed.dat", "w") as file:
        file.writelines(header)
        for x, y in data:
            file.write(f"{x}\t{y}\n")

# def add_flap_effect(parameter):
#     # airfoil.datを読み込み
#     with open("C:/CFD2/airfoil.dat", "r") as file:
#         lines = file.readlines()

#     # ヘッダーを保持
#     header = lines[:2]
#     data_lines = lines[2:]

#     # 座標データを取得
#     x = []
#     y = []
#     for line in data_lines:
#         x_str, y_str = line.strip().split()
#         x_val = float(x_str)
#         y_val = float(y_str)
#         x.append(x_val)
#         y.append(y_val)

#     # x >= 0.8 のデータを抽出
#     flap_indices = [i for i, val in enumerate(x) if val >= 0.8]
#     flap_x = [x[i] for i in flap_indices]
#     flap_y = [y[i] for i in flap_indices]

#     # x = 0.8 に最も近い、かつ y < 0 の点を回転中心とする
#     valid_indices = [i for i in flap_indices if y[i] < 0]
#     if valid_indices:
#         rotation_index = min(valid_indices, key=lambda i: abs(x[i] - 0.8))
#         cx, cy = x[rotation_index], y[rotation_index]
#     else:
#         raise ValueError("No valid points found for rotation center.")

#     # 回転角度（時計回りにparameter度）
#     theta = np.radians(-parameter)

#     # 回転行列を適用する関数
#     def rotate_point(x, y, cx, cy, theta):
#         x_shifted = x - cx
#         y_shifted = y - cy
#         x_rot = x_shifted * np.cos(theta) - y_shifted * np.sin(theta) + cx
#         y_rot = x_shifted * np.sin(theta) + y_shifted * np.cos(theta) + cy
#         return x_rot, y_rot

#     # フラップデータに回転を適用
#     for i in flap_indices:
#         x[i], y[i] = rotate_point(x[i], y[i], cx, cy, theta)

#         # 翼全体を回転・拡大して x が最大の点を (1, 0) に移動
#     max_index = np.argmax(x)
#     max_x, max_y = x[max_index], y[max_index]

#     # 翼前縁（x = 0, y = 0）を回転中心とする
#     cx, cy = 0.0, 0.0

#     # 移動のための回転角度を計算
#     theta_adjust = np.arctan2(-max_y, max_x)

#     # 拡大係数を計算
#     scale_factor = 1 / max_x

#     # 回転と拡大を適用する関数
#     def rotate_and_scale(x, y, cx, cy, theta, scale):
#         x_shifted = x - cx
#         y_shifted = y - cy
#         x_rot = (x_shifted * np.cos(theta) - y_shifted * np.sin(theta)) * scale + cx
#         y_rot = (x_shifted * np.sin(theta) + y_shifted * np.cos(theta)) * scale + cy
#         return x_rot, y_rot

#     # 翼全体に回転と拡大を適用
#     for i in range(len(x)):
#         x[i], y[i] = rotate_and_scale(x[i], y[i], cx, cy, theta_adjust, scale_factor)
    
#     x[-1]=1
#     y[-1]=0

#     # データを書き込み
#     with open("C:/CFD2/work/airfoil-snowed.dat", "w") as file:
#         file.writelines(header)
#         for xi, yi in zip(x, y):
#             file.write(f"{xi}\t{yi}\n")
    
#     return theta_adjust

def compile_mesh():
    command = ["gfortran", "-o", "build/tmesh", "C:/CFD2/mesh/t_mesh.f", "C:/CFD2/mesh/mesh_io.f", "C:/CFD2/mesh/mesh_cv.f"]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode == 0:
        print("Compilation successful")
        execute_tmesh()
    else:
        print("Compilation failed")
        print(result.stderr)

def execute_tmesh():
    result = subprocess.run(["./build/tmesh"], capture_output=True, text=True)
    if result.returncode == 0:
        print("Execution successful")
        print(result.stdout)
        compile_euler()
    else:
        print("Execution failed")
        print(result.stderr)

def compile_euler():
    command = ["gfortran", "-o", "build/2deuler", "C:/CFD2/euler/main.f","C:/CFD2/euler/bc.f","C:/CFD2/euler/io.f"]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode == 0:
        print("Euler compilation successful")
    else:
        print("Euler compilation failed")
        print(result.stderr)

def create_mesh_d(fsmach, alpha):
    with open("C:/CFD2/work/mesh.d", "w") as file:
        file.write("&INPUT\n")
        file.write(f"   fsmach={fsmach}, alpha={alpha},\n")
        file.write("   mesh_file='work/mesh.xyz', q_file='work/mesh.q',\n")
        file.write("   restart=.fault.,\n")
        file.write("   n_steps=20000, cdt=0.9, eps4=20.,  jtail=17\n")
        file.write(" &END\n")



def execute_euler():
    process = subprocess.Popen(["./build/2deuler"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = [], []

    for line in iter(process.stdout.readline, ''):
        print(line, end='')
        stdout.append(line)
    process.stdout.close()

    for line in iter(process.stderr.readline, ''):
        print(line, end='')
        stderr.append(line)
    process.stderr.close()

    process.wait()

    if process.returncode == 0:
        print("Euler execution successful")
        
        # 出力を解析
        output = ''.join(stdout)
        cl_integration = re.search(r"Cl \(integration of Cp_lower - Cp_upper\) =\s+([\d.E+-]+)", output)
        cl_surface = re.search(r"Cl \(surface integration of pressure\) =\s+([\d.E+-]+)", output)
        cd_surface = re.search(r"Cd \(surface integration of pressure\) =\s+([\d.E+-]+)", output)
        
        results = {
            "Cl_integration": float(cl_integration.group(1)) if cl_integration else None,
            "Cl_surface": float(cl_surface.group(1)) if cl_surface else None,
            "Cd_surface": float(cd_surface.group(1)) if cd_surface else None
        }
        return results
    else:
        print("Euler execution failed")
        print(''.join(stderr))

def save_result(parameter, fsmach, alpha, results, output_dir, append=False):
    results.update({"parameter": parameter, "fsmach": fsmach, "alpha": alpha})
    output_file = os.path.join(output_dir, "result.csv")
    mode = 'a' if append else 'w'
    
    # CSVファイルに結果を書き込む
    with open(output_file, mode, newline='') as csvfile:
        fieldnames = ["parameter", "fsmach", "alpha", "Cl_integration", "Cl_surface", "Cd_surface"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        if not append:
            writer.writeheader()
        writer.writerow(results)
    
    # フォルダを作成してファイルをコピー
    folder_name = os.path.join(output_dir, f"result_{fsmach}_{alpha}")
    os.makedirs(folder_name, exist_ok=True)
    shutil.copy("work/mesh.q", folder_name)
    shutil.copy("work/residual.dat", folder_name)
    shutil.copy("work/surf.dat", folder_name)

if __name__ == "__main__":
    fsmach_values = [0.7, 0.8, 0.9]
    alpha_values = [0.0,5.0, 10.0]
    parameters = [0.0, 0.5, 1]
    # ここからメイン処理

    create_directories()

    first_write = True

    # 出力フォルダのベースパスを設定
    base_output_dir = "output"


    for parameter in parameters:
        # フラップ角に対応するフォルダを作成
        parameter_folder = os.path.join(base_output_dir, f"flap_angle_{parameter}")
        os.makedirs(parameter_folder, exist_ok=True)

        add_snow_effect(parameter=parameter)
        # theta_adjust = add_flap_effect(parameter=parameter)
        compile_mesh()
        print(f"積雪パラメータ：{parameter}")

        for fsmach in fsmach_values:
            for alpha in alpha_values:
                create_mesh_d(fsmach=fsmach, alpha=alpha) #+ theta_adjust)
                results = execute_euler()

                # save_result関数に出力ディレクトリを渡す
                save_result(
                    parameter=parameter,
                    fsmach=fsmach,
                    alpha=alpha,
                    results=results,
                    output_dir=parameter_folder,
                    append=not first_write
                )
                first_write = False

    os.remove("work/mesh.q")
    os.remove("work/residual.dat")
    os.remove("work/surf.dat")
    exit(0)
