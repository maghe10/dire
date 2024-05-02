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

#move to storage
storage = '/root/sequencing/storage/unicycler/conservative/sample' + sample
output = '/root/sequencing/intermediate/unicycler/conservative_to_normal/sample' + sample

command = "mkdir -p " + output
print(command)
os.system(command)

command = "cp " + storage + "/* " + output
print(command)
os.system(command)
  