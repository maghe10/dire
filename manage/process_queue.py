import os
import socket

hostname = socket.gethostname()
manage_root = "/root/dire/data/manage/" + hostname

queuedir = manage_root + "/queue/"
queueddir = manage_root + "/queued/"
finisheddir = manage_root + "/finished/"

os.system("mkdir -p " + queuedir)
os.system("mkdir -p " + queueddir) 
os.system("mkdir -p " + finisheddir)

print("Processing queue in "+ queuedir)
listing = sorted(os.listdir(queuedir))
print("Number of scripts in queue: " + str(len(listing)))
if (len(listing)>0):
	print(listing[0])

while len(listing)>0:
	#pick first in queue
	first = listing[0]
	#move to queued
	os.rename(queuedir + "/" + first, queueddir + "/" + first)
	#execute first item
	rv = os.system("python " + queueddir + "/" + first)
	if rv != 0:
          break
	#move to finished
	os.rename(queueddir + "/" + first, finisheddir + "/" + first)
	#relist
	listing = sorted(os.listdir(queuedir))


      
