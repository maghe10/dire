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

#mkdir -p /root/sequencing/intermediate/quest/sample2
#quast.py /root/sequencing/intermediate/spades/sample2/contigs.fasta -o /root/sequencing/intermediate/quast/sample2
contigs = '/root/sequencing/storage/assembly/spades_optimal/sample' + sample + ".fasta"
path_to_program = 'quast.py'
output = '/root/sequencing/intermediate/quast/spades_optimal/sample' + sample
os.makedirs("output", exist_ok=True)

invoke = path_to_program + " --gene-finding -o " + output + " " + contigs
print(invoke)
result = "skip" 
result = subprocess.run([invoke],shell=True)
print(result),
