import subprocess
import argparse
import os

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-s", "--sample", default=None)
args = parser.parse_args()

if args.sample == None:
  print("Usage: -s <sample>")
  raise SystemExit

#################
sample = args.sample
print("\nSample: ")
print(sample)

aReadForward = "sample" + sample + "_1_val_1.fq.gz"
aReadBackward = "sample" + sample + "_2_val_2.fq.gz"
output = '/root/sequencing/intermediate/fastqtrimmed'

os.system("mkdir -p " + output) 

## Fastqc 
f = '/root/sequencing/storage/trimmed/' + aReadForward
b = '/root/sequencing/storage/trimmed/' + aReadBackward
path_to_program = 'fastqc'
args =  "-o" +  " " + output + " " + f + " " + b
invoke = path_to_program + " " + args
result = "skip" 
result = subprocess.run([invoke],shell=True)
print(result)
  