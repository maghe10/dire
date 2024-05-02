import subprocess
import argparse

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

aReadForward = "sample" + sample + "_1.fastq.gz"
aReadBackward = "sample" + sample + "_2.fastq.gz"
output = '/root/sequencing/intermediate/fastq'

## Fastqc 
f = '/root/sequencing/in/reads/' + aReadForward
b = '/root/sequencing/in/reads/' + aReadBackward
path_to_program = '/usr/bin/fastqc'
args =  "-o" +  " " + output + " " + f + " " + b
invoke = path_to_program + " " + args
result = "skip" 
result = subprocess.run([invoke],shell=True)
print(result)
  