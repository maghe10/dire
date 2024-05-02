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


with open("trimgalore_dire.py") as f:
    exec(f.read())

with open("spades_dire.py") as f:
    exec(f.read())

with open("filter_dire.py") as f:
    exec(f.read())

with open("confindr_dire.py") as f:
    exec(f.read())

with open("quast_dire.py") as f:
    exec(f.read())

with open("store_dire.py") as f:
    exec(f.read())
