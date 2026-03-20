import sys
import torch

print("=== ENVIRONMENT TEST ===")

print("Python version:", sys.version)
print("PyTorch version:", torch.__version__)
print("CUDA available:", torch.cuda.is_available())
print("CUDA device count:", torch.cuda.device_count())

# CPU test
try:
    x = torch.randn(2, 2)
    y = torch.randn(2, 2)
    print("\nCPU test OK:", x @ y)
except Exception as e:
    print("CPU test FAILED:", e)

# GPU test
if torch.cuda.is_available():
    try:
        device = torch.device("cuda")
        x = torch.randn(2, 2).to(device)
        y = torch.randn(2, 2).to(device)
        print("\nGPU device:", torch.cuda.get_device_name(0))
        print("GPU test OK:", x @ y)
    except Exception as e:
        print("GPU test FAILED:", e)
else:
    print("\nGPU not available → GPU test skipped.")
