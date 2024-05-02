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

#filter optimal spades assembly, At least 500 bp
assembly = '/root/sequencing/storage/unicycler/conservative/sample' + sample + "/assembly.fasta"
path_to_program = 'seqkit'
outputdir = '/root/sequencing/storage/assembly_unicycler_conservative'
output = outputdir + "/sample"+ sample + ".fasta"

invoke = "mkdir -p " + outputdir
print(invoke)
result = "skip" 
result = subprocess.run([invoke],shell=True)
print(result),

invoke = path_to_program + " seq -m 500 " + assembly + " > " + output
print(invoke)
result = "skip" 
result = subprocess.run([invoke],shell=True)
print(result),
