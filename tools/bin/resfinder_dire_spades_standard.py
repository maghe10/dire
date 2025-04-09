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


input = '/root/sequencing/storage/assembly/spades_standard/sample' + sample + '.fasta'
print(input)

outputdir = '/root/sequencing/storage/resfinder/spades_standard'
output = outputdir + "/sample"+ sample + "/"

invoke = "mkdir -p " + output
print(invoke)
result = "skip" 
result = subprocess.run([invoke],shell=True)
print(result),


invoke = "python3  -m resfinder -o " + output + " -s \"Escherichia coli\" -l 0.6 -t 0.8 --acquired --point -db_res /root/resfinder/db_resfinder -ifa " + input 
print(invoke)
result = "skip" 
result = subprocess.run([invoke],shell=True)
print(result),

