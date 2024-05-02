import subprocess
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-s", "--sample", default=None)
args = parser.parse_args()

if args.sample == None:
  print("Usage: -s <sample>")
  raise SystemExit

## assemble with spades 
## python /root/sequencing/spades/SPAdes-3.15.4/bin/spades.py -k 21,33,55,77 --isolate
# SPAdes with --isolate
sample = args.sample
print("\nSample: ")
print(sample)

f = '/root/sequencing/intermediate/trimmed/sample' + sample + "/sample" + sample + "_1_val_1.fq.gz"
b = '/root/sequencing/intermediate/trimmed/sample' + sample + "/sample" + sample + "_2_val_2.fq.gz"
path_to_program = '/root/sequencing/spades/SPAdes-3.15.4/bin/spades.py'
output = '/root/sequencing/intermediate/spades/sample' + sample
invoke = 'python' + " " + path_to_program + " -m 10 -t 4 -k 21,33,55,77 --isolate --cov-cutoff auto " + "-1 " + f + " -2 " + b + " -o " + output
print(invoke)
result = "skip" 
result = subprocess.run([invoke],shell=True)
print(result)
