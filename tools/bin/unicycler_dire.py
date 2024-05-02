import subprocess
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-s", "--sample", default=None)
args = parser.parse_args()

if args.sample == None:
  print("Usage: -s <sample>")
  raise SystemExit

sample = args.sample
print("\nSample: ")
print(sample)

trimmedRoot = "/root/sequencing/storage/trimmed"

f = trimmedRoot + "/sample" + sample + "_1_val_1.fq.gz"
b = trimmedRoot + "/sample" + sample + "_2_val_2.fq.gz"
output = '/root/sequencing/intermediate/unicycler/conservative/sample' + sample

path_to_program = '/root/miniconda3/envs/dire/bin/unicycler'
invoke = 'python' + " " + path_to_program  + " -1 \'" + f + "\' -2 \'" + b + "\' -o " + output + " --mode conservative"
print(invoke)
result = "skip" 
result = subprocess.run([invoke],shell=True)
print(result)
