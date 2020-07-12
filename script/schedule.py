#!/usr/bin/env python
import subprocess
import time
import sys
import os

# a = subprocess.Popen("ssh bsidna2 /mnt/cluster/xdeng/intense.py 100000000", shell=True)
# b = subprocess.Popen("ssh bsidna3 /mnt/cluster/xdeng/intense.py 100000000", shell=True)

# while a.poll() is None or b.poll() is None:
    # time.sleep(3)
    # print 'a', a.pid, a.poll(),'b', b.pid, b.poll()

def job_done(job):
	parts = job.strip().split()
	filename=''
	for part in parts:
		if 'xml' in part:
			filename= part
			break
	if os.path.isfile(filename) and os.stat(filename).st_size > 0:
		return True
	else:
		return False

def readdata(file, dic=False):
	if dic==False: jobs=[]
	else: jobs={}
	f=open(file, 'r')
	for line in f:
		if dic==False:
			jobs.append(line.strip())
		else:
			jobs[line.strip()]=None
	f.close()
	return jobs

def run(i):
	j=0
	running=0
	while j < len(servernames) and i< len(jobs):
		job=jobs[i]
		s=servernames[j]
		ps=servers[s]
		status=''
		if ps is not None:
			status = ps.poll()
		if not job_done(job) and (ps is None or status is not None): #server and job ready
			ps = subprocess.Popen('ssh '+s+' '+job, shell=True)
			servers[s]=ps
			i+=1 #next job
			j+=1 #next server
		elif not job_done(job): #job ready server not ready
			j+=1 #next server
		elif ps is None or status is not None: #job not ready server ready
			i+=1 #next job
		else: #job server both not ready
			i+=1 #next job
			j+=1 #next server
	for (server,ps) in servers.items():
		if ps is not None and ps.poll() is None:
			running+=1
	return i, running

def run_all():
	i, running= run(0) #start from job 0 on all servers
	while running >0:
		print 'i', i, 'of', len(jobs), 'jobs. running servers', running
		time.sleep(15)
		i, running = run(i)

if __name__ == "__main__":
	jobs = readdata(sys.argv[1])
	servers=readdata(sys.argv[2], True)
	servernames=servers.keys()
	run_all()
