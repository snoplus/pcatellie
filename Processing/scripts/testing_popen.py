import subprocess
import shlex
import time

cmd = "du -sh /home/michal/"
args = shlex.split(cmd)

process = subprocess.Popen( args, stdout=subprocess.PIPE, stderr=subprocess.PIPE )
#out, err = process.communicate()

while True:
	#print out
	print process.poll()
	print " "
	time.sleep(0.1)
