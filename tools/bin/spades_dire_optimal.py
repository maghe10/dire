import subprocess
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-s", "--sample", default=None)
args = parser.parse_args()

if args.sample == None:
  print("Usage: -s <sample>")
  raise SystemExit

## assemble with spades 
## python /root/sequencing/spades/SPAdes-3.15.4/bin/spades.py -k 21,33,55,77 --isolate
# SPAdes with --isolate
sample = args.sample
print("\nSample: ")
print(sample)

# Run to find optimal kmere length 
# find /root/sequencing/storage/unicycler/ -type f -name "unicycler.log" -exec grep -Hn "best$" '{}' ';' | grep -v 127

kmerlength = "27,53,71,87,99,111,119,127"
if sample  == "19":
  kmerlength = "27,53,71,87,99"
if sample  == "25":
  kmerlength = "27,53,71,87,99,111"
if sample  == "38":
  kmerlength = "27,53"
if sample  == "46":
  kmerlength = "27,53,71,87"
if sample  == "52":
  kmerlength = "27,53,71,87,99,111,119"
if sample  == "54":
  kmerlength = "27,53,71,87,99,111"
if sample  == "56":
  kmerlength = "27,53,71,87,99"
if sample  == "58":
  kmerlength = "27,53,71,87,99,111,119"
if sample  == "61":
  kmerlength = "27,53,71,87"
if sample  == "76":
  kmerlength = "27,53,71,87,99,111,119"
if sample  == "97":
  kmerlength = "27,53,71,87,99,111,119"
if sample  == "116":
  kmerlength = "27,53,71,87,99,111"

print("\nkmer length: ")
print(kmerlength)


f = "/root/sequencing/storage/trimmed/sample" + sample + "_1_val_1.fq.gz"
b = "/root/sequencing/storage/trimmed/sample" + sample + "_2_val_2.fq.gz"
path_to_program = '/root/sequencing/spades/SPAdes-3.15.4/bin/spades.py'
output = '/root/sequencing/intermediate/spadesoptimal/sample' + sample
invoke = 'python' + " " + path_to_program + " -m 10 -t 4 -k " + kmerlength +" --isolate --cov-cutoff auto " + "-1 " + f + " -2 " + b + " -o " + output
print(invoke)
result = "skip" 
result = subprocess.run([invoke],shell=True)
print(result)
