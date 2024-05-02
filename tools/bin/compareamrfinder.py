import os
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

storagedir = "/root/sequencing/storage/"
tmpdir = "/root/sequencing/tmp/"
diffdir = tmpdir + "amrfinder/"


os.system("mkdir -p " + diffdir)

for sample in samples:
	#awk '{ print $6 }' amrfinder/sample2.tsv > amrfinderoptimal_2_genes.txt
	command="awk \'{ print $6 }\' " +storagedir + "amrfinder/sample" + sample + ".tsv > "+ diffdir + "amrfinder_" + sample + "_genes.txt"
	#print(command)
	os.system(command)
	command="awk \'{ print $6 }\' " +storagedir + "amrfinder_optimal/sample" + sample + ".tsv > "+ diffdir + "amrfinderoptimal_" + sample + "_genes.txt"
	#print(command)
	os.system(command)
	command = "sort " + diffdir + "amrfinder_" + sample + "_genes.txt > " + diffdir + "amrfinder_" + sample + "_genes_sortes.txt"
	#print(command)
	os.system(command)
	command = "sort " + diffdir + "amrfinderoptimal_" + sample + "_genes.txt > " + diffdir + "amrfinderoptimal_" + sample + "_genes_sortes.txt"
	#print(command)
	os.system(command)
	command = "diff " +  diffdir + "amrfinderoptimal_" + sample + "_genes_sortes.txt" + " " + diffdir + "amrfinder_" + sample + "_genes_sortes.txt" + ">" + diffdir + "diff_" + sample + "armfinder.txt"
	#print(command)
	os.system(command)
	command = "cat " + diffdir + "diff_" + sample + "armfinder.txt"
	os.system(command)

