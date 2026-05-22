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

aRead = "sample" + sample + ".fastq.gz"
output = '/root/sequencing/intermediate/trimmed/sample' + sample

## Trim with trim galore
r = '/root/dire/data/Illumina/singleend/' + aRead
path_to_program = '/root/sequencing/trimgalore/TrimGalore-0.6.10/trim_galore'
args = r + " -o" +  " " + output
invoke = path_to_program + " " + args
result = "skip" 
result = subprocess.run([invoke],shell=True)
print(result)
  