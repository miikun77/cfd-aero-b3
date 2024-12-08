import subprocess
import os
import re
import csv

def create_directories():
    os.makedirs("build", exist_ok=True)
    os.makedirs("work", exist_ok=True)

def add_snow_effect(parameter):
    # airfoil.datを読み込み
    with open("airfoil.dat", "r") as file:
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
    with open("work/airfoil-snowed.dat", "w") as file:
        file.writelines(header)
        for x, y in data:
            file.write(f"{x}\t{y}\n")

def compile_mesh():
    command = ["gfortran", "-o", "build/tmesh", "mesh/t_mesh.f", "mesh/mesh_io.f", "mesh/mesh_cv.f"]
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
    command = ["gfortran", "-o", "build/2deuler", "euler/main.f","euler/bc.f","euler/io.f"]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode == 0:
        print("Euler compilation successful")
    else:
        print("Euler compilation failed")
        print(result.stderr)

def create_mesh_d(fsmach, alpha):
    with open("work/mesh.d", "w") as file:
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
        
        # CSVに保存
        with open("result.csv", "w", newline='') as csvfile:
            fieldnames = ["Cl_integration", "Cl_surface", "Cd_surface"]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            
            writer.writeheader()
            writer.writerow(results)
    else:
        print("Euler execution failed")
        print(''.join(stderr))

if __name__ == "__main__":
    create_directories()
    create_mesh_d(fsmach=0.800, alpha=1.25)
    add_snow_effect(parameter=0.05)
    compile_mesh()
    execute_euler()
