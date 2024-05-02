import subprocess
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--infile", default=None)
parser.add_argument("-o", "--outfile", default=None)
args = parser.parse_args()

if args.outfile == None:
  print("Usage: -i <in> -o <out>")
  raise SystemExit

#if args.infile == None:
#  print("Usage: -i <in> -o <out>")
#  raise SystemExit

print("\nInfile " + args.infile + "\nOutfile " + args.outfile + "\n")



### awk '/^S/{print ">"$2"\n"$3}' in.gfa | fold > out.fa
command = "awk \'/^S/{print \">\"$2\"\\n\"$3}\' " + args.infile + " | fold > " + args.outfile
print("command: " + command)
 
result = subprocess.run([command],shell=True)
print(result),
