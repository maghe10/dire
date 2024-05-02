import os
import subprocess
import argparse

samplesA = ["2","3","4","5","6","8"]
samplesB = ["9","11","12","13","14","18","19","22"] 
samplesC = ["24","25","26","28","29","34","35","36"]
samplesD = ["39","40","43","44","45","47","48","50"]
samplesE = ["51","52","56","57","59","61","62","64"]
samplesF = ["66","67","69","70","72","76","80","82"]
samplesG = ["83","84","86","87","88","89","91","93"]
samplesH = ["94","96","97","98","101","102","103","105"] 
samplesI = ["109","111","112","113","120"] 
samplesJ = ["17","20","23","27","32","37","38","41","46"]  
samplesK = ["49","53","54","58","75","95","100","115","116"]
samplesL = ["21","63","77","78","114"]
samplesM = ["121","122","123","124","55B","60B","65B","79B","85B","110B"]
samplesN = ["125"]
allSamples = samplesA + samplesB + samplesC + samplesD + samplesE + samplesF + samplesG + samplesH + samplesI + samplesJ + samplesK + samplesL + samplesM + samplesN


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-s", "--sample", default=None)
parser.add_argument("-b", "--batch", default=None)
args = parser.parse_args()

samples = allSamples
if args.sample != None:
	sample = args.sample
	print("\nSample: ")
	print(sample)
	samples = [sample]

if args.batch !=None:
	batch = args.batch
	print("\nBatch: ")
	print(batch)
	samples = globals()['samples'+batch]

#################

print(samples)



os.chdir("/root/sequencing/tools/bin")

tmpdir = "/root/sequencing/storage/tmp/"

quastdir = "/root/sequencing/storage/old_withspades/"

for sample in samples:
	finddir = quastdir + "sample" + sample + "/" 
	targetdir = tmpdir + "sample" + sample + "/"
	os.system("mkdir -p " + targetdir)
	command = "find " + finddir + " -type f -name report.tsv -exec cp '{}' "+targetdir +" ';'"
	print(command)
	result = os.system(command)
	print(result)	
	command = "find " + finddir + " -type f -name transposed_report.tsv -exec cp '{}' "+targetdir +" ';'"
	print(command)
	result = os.system(command)
	print(result)	


fastqcdir= "/root/sequencing/storage/trimmed_fastqc/"

for sample in samples:
	finddir = fastqcdir + "sample" + sample + "/"
	targetdir = tmpdir + "sample" + sample + "/"
	command = "find " + finddir + " -type f \! -name \*gz -exec cp '{}' "+targetdir +" ';'"
	print(command)
	result = os.system(command)
	print(result)	

