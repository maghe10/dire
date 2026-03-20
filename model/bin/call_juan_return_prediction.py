from pathlib import Path
import subprocess

# Full path to Rscript.exe
rscript = r"C:\Users\xhessm\AppData\Local\Programs\R\R-4.5.2\bin\Rscript.exe"

# Project root (working directory)
project_root = Path("..") / "Confidence-based-Prediction-of-Antibiotic-Resistance"

# Path to the R script
script = project_root / "CICP" / "return_prediction_sets.R"

# Arguments you want to pass to the R script
args = ["data/input_with_patient_data.csv", "data/input_without_patient_data.csv", "full"]   # <-- replace with your real args

# Build the command
cmd = [rscript, str(script)] + args

# Run the script
result = subprocess.run(
    cmd,
    cwd=str(project_root),     # R script runs here
    capture_output=True,
    text=True
)

print("STDOUT:\n", result.stdout)
print("STDERR:\n", result.stderr)
print("Exit code:", result.returncode)
