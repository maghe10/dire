import os


queuedir = "/root/sequencing/tools/manage/gu/queue/"
queueddir = "/root/sequencing/tools/manage/gu/queued/"
finisheddir = "/root/sequencing/tools/manage/gu/finished/"

listing = sorted(os.listdir(queuedir))
print(len(listing))
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


      
