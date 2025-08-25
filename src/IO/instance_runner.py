import os
import re
import subprocess
import time

class InstanceRunner:
    root_path: str
    instance_path: str

    def __init__(self, parameters:dict, root_path: str, instance_path: str, header_path: str):
        self.parameters = parameters
        self.root_path = root_path
        self.instance_path = instance_path
        self.header_path = header_path

    def compile(self):
        try:
            subprocess.run(
                ["make", "prd"],
                cwd=self.root_path,
                check=True
            )
            print("Compilation finished successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Compilation failed with error code {e.returncode}")
        except FileNotFoundError:
            print("`make` command not found. Please ensure it is installed and in PATH.")

    def get_instance_columns(self):
        with open(self.instance_path) as f:
            return int(f.readline().split()[1])
        
    def get_header_columns(self):
        with open(self.header_path) as f:
            for line in f:
                match = re.match(r"#define\s+COLUMNS\s+(\d+)", line)
                if match:
                    return int(match.group(1))
        raise ValueError("No '#define COLUMNS <number>' found in instance file.")
    
    def replace_headers_columns(self, columns: int):
        with open(self.header_path, "r") as f:
            header = f.read()

        header = re.sub(r'\#define\s*COLUMNS\s*\d+', f"#define COLUMNS {columns}", header)

        with open(self.header_path, 'w') as f:
            f.write(header)

    def call_solver(self):
        parameters = [f"--{key}={val}" for key, val in self.parameters.items()]
        command = ["./main_prd", f"--filePath={self.instance_path}"] + parameters
        start_time = time.time()
        result = subprocess.run(
            command,
            cwd=self.root_path,
            check=True,
            capture_output=True,
            text=True
        )
        end_time = time.time()
        return result.stdout, end_time - start_time

    def run(self):
        instance_columns = self.get_instance_columns()
        header_columns = self.get_header_columns()
        if not instance_columns == header_columns:
            self.replace_headers_columns(instance_columns)
            self.compile()
        return self.call_solver()
