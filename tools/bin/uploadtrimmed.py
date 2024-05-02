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
output = '/root/sequencing/intermediate/trimmed/sample' + sample
storage = '/root/sequencing/storage/trimmed'
command = "mv " + output + "/* " + storage
print(command)
os.system(command)

command = "rmdir " + output
print(command)
os.system(command)
