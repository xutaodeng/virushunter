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
	global jobStatus
	parts = job.strip().split()
	#filename=parts[-1]+'.daa' #diamond
	for part in parts:
		if '.xml' in part or part.endswith('.m8') or part.endswith('.sam'):
			filename= part
			break
	if jobStatus[job]is not None and (os.path.isfile(filename) and os.stat(filename).st_size > 0) or (os.path.isfile(filename+'.xml') and os.stat(filename+'.xml').st_size > 0):
		return True
	else:
		return False

def readdata(file):
	jobs=[]
	f=open(file, 'r')
	for line in f:
		jobs.append(line.strip())
	f.close()
	return jobs

def initiate():
	global jobStatus
	for job in jobs:
		jobStatus[job]=None

def run(i):
	global badjobs
	global jobStatus
	j=0
	running=set([])
	donejobs=0
	
	for (server,info) in servers.items():
		ps, lastjob=info
		if ps is not None and ps.poll() is None:
			running.add(lastjob)
			
	while j < len(servernames) and i< len(jobs):
		job=jobs[i]
		s=servernames[j]
		ps, lastjob=servers[s]
		jobname=jobs[lastjob]
		if ps is not None:
			jobStatus[jobname]=ps.poll()

		physicalserver=s.rsplit('_',1)[0]
		status=''
		
		if job_done(job) or (i in running) or (i in badjobs): #next job
			i+=1
			continue

		if ps is not None: #check server
			status = ps.poll()
			if status is not None and status!=0 and not job_done(jobs[lastjob]): #last job finished and bad
				badjobs.append(lastjob)
				badstatus.append(status)
		if  ps is None or status is not None: #server ready
			# if status is not None:
				# print status, lastjob
			ps = subprocess.Popen('ssh '+physicalserver+' '+job, shell=True)
			running.add(i)
			servers[s]=(ps, i)
			i+=1 #next job
			j+=1 #next server
		else: 
			j+=1 #next server

	ind, badind=0, 0
	for kkk in badjobs:
		if not job_done(jobs[kkk]):
			print 'status=',badstatus[ind], 'jobno=', kkk, jobs[kkk]
			badind+=1
		ind+=1
	if badind==0: badjobs=[]
		
	for job in jobs:
		if job_done(job):
			donejobs+=1

	return i, len(running), donejobs, len(badjobs)

def run_all():
	i=0
	i, runy, donejobs, bad = run(i) #start from job 0 on all servers
	while runy >0 and donejobs<len(jobs):
		print i, 'of', len(jobs), 'jobs. running servers=', runy, 'donejobs=', donejobs, 'badjobs=', bad, jobs[0].split()[0]
		time.sleep(20)
		i, runy, donejobs, bad= run(i)
	print i, 'of', len(jobs), 'jobs. running servers=', runy, 'donejobs=', donejobs, 'badjobs=', bad, jobs[0].split()[0]
	time.sleep(20)
	

if __name__ == "__main__":
	print 'running...'
	badjobs=[]
	badstatus=[]
	jobs = readdata(sys.argv[1])
	serverslist=readdata(sys.argv[2])
	jobStatus={}
	try: serverload=int(sys.argv[3])
	except: serverload=1
	servers={}
	for server in serverslist:
		for i in xrange(serverload):
			servers[server+'_'+str(i)]=(None, -1)
	servernames=servers.keys()
	initiate()
	run_all()
	run_all() #make up the failed runs
	
