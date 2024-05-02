from datetime import datetime
import subprocess
import argparse
import os

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-s", "--sample", default=None)
args = parser.parse_args()

if args.sample == None:
  print("Usage: -s <sample>")
  raise SystemExit

sample = args.sample
print("\nSample: ")
print(sample)

intermediatedir = '/root/sequencing/intermediate'

storagedir = '/root/sequencing/storage/quast/sample' + sample

result = os.system("mkdir -p " + storagedir)

command = "mv " +  intermediatedir + "/quast/sample" + sample + "/* "+ storagedir
print(command)
os.system(command)

command = "rmdir " + intermediatedir + "/quast/sample" + sample
print(command)
os.system(command)
