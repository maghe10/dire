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

# Find read files.
# Run confindr. This assumes that you have already downloaded the databases. If you haven't,
# you can run confindr.check_for_databases_and_download(database_location='path/where/you/want/to/download, tmpdir='a/tmp/dir')
pair = ['/root/sequencing/storage/trimmed/sample' + sample + '_1_val_1.fq.gz','/root/sequencing/storage/trimmed/sample' +  sample + '_2_val_2.fq.gz']
print(pair)
confindr.find_contamination(pair=pair,
                            	forward_id='_1', # change if yours is different
                            	threads=4,
				base_cutoff=3,	 
                                output_folder='/root/sequencing/intermediate/confindr/sample' + sample,
                                databases_folder='/root/.confindr_db')