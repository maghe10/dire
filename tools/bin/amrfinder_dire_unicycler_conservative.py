from confindr_src import confindr
#import confindr
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


input = '/root/sequencing/storage/assembly/unicycler_conservative/sample' + sample + '.fasta'
print(input)

outputdir = '/root/sequencing/storage/amrfinder/unicycler_conservative'
output = outputdir + "/sample"+ sample + ".tsv"

invoke = "mkdir -p " + outputdir
print(invoke)
result = "skip" 
result = subprocess.run([invoke],shell=True)
print(result),


invoke = "amrfinder --plus -n " + input + " -O Escherichia > " + output 
print(invoke)
result = "skip" 
result = subprocess.run([invoke],shell=True)
print(result),
