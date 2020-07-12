#!/usr/bin/env python
import os
import sys

f=open(sys.argv[1], 'r') #input script
of=open(sys.argv[2],'w') #output parallel script

cwd = os.getcwd()

servers=[ 'bsidna2', 'bsidna3', 'bsidna4', \
 'bsidna5', 'bsidna6','bsidna7','bsidna8', \
  'bsidna9', 'bsidna10','bsidna11', 'bsidna12','bsidna13', \
  'bsidna14', 'bsidna15', 'bsidna16', 'bsidna17', \
'bsidna18', 'bsidna19', 'bsidna20', \
 'bsidna22', 'bsidna23', 'bsidna24', 'bsidna25', 'bsidna26', \
'bsidna27','bsidna28','bsidna29'
]
nservers=len(servers)
job=0

for line in f:
	if line.strip().startswith('#'): continue
	if line.strip()=='': continue
	serverTag = 'ssh '+servers[job%nservers]
	if '"' in line and "'" in line:
		print 'line contains both "'+" and ', fix it!!!" 
	if '"' in line:
		of.write(serverTag+" 'cd "+cwd+' && '+ line.strip()+" '&\n")
	else:
		of.write(serverTag+' "cd '+cwd+' && '+ line.strip()+' "&\n')
	job+=1
	if job%nservers==0:  of.write('wait\n')
		
f.close()
of.close()