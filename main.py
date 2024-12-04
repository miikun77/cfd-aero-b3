import subprocess
import os

def create_directories():
    os.makedirs("build", exist_ok=True)
    os.makedirs("work", exist_ok=True)

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
        execute_euler()
    else:
        print("Euler compilation failed")
        print(result.stderr)

def execute_euler():
    result = subprocess.run(["./build/2deuler"], text=True)
    if result.returncode == 0:
        print("Euler execution successful")
    else:
        print("Euler execution failed")

# ...existing code...

if __name__ == "__main__":
    create_directories()
    compile_mesh()
    # ...existing code...
