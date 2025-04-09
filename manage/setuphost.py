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


      
