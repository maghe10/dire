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

#filter spades assembly, At least 500 bp
scaffolds = '/root/sequencing/intermediate/spades/sample' + sample + "/scaffolds.fasta"
path_to_program = 'seqkit'
outputdir = '/root/sequencing/intermediate/assembly/sample' + sample
output = outputdir + "/sample"+ sample + ".fasta"

invoke = "mkdir -p " + outputdir
print(invoke)
result = "skip" 
result = subprocess.run([invoke],shell=True)
print(result),

invoke = path_to_program + " seq -m 500 " + scaffolds + " > " + output
print(invoke)
result = "skip" 
result = subprocess.run([invoke],shell=True)
print(result),
