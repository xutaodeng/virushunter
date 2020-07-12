#!/usr/bin/env python
import sys
import fcntl

ff=open(sys.argv[1], 'r')
f=open(sys.argv[2], 'r')
of=open(sys.argv[3],'w') #result
of2=open(sys.argv[4],'a') #log file

cluster={}
for line in ff:
	if line.startswith('CLUSTER'):
		cnumber=line.strip().split()[1]
	elif line.startswith('FAMILIES'):
		family=' '.join(line.strip().split()[1:])
	elif line.startswith('GENERA'):
		genera=' '.join(line.strip().split()[1:])
		cluster[cnumber]=genera+'&'+family
ff.close()

hits=0
for line in f:
	if line.startswith('#'):
		if line.startswith('#='):
			of.write(line+'\n')
			if line.startswith('#=GS'):
				hits+=1
		else:
			continue
	elif 'vFam_' in line:
		parts=line.strip().split()
		contig=parts[0]
		cnumber=parts[2].split('_')[1]
		EVAL=parts[4]
		annot=cluster[cnumber]
		out=contig+' Evalue:'+EVAL+' '+annot
		of.write(out+'\n')

f.close()
of.close()
fcntl.flock(of2, fcntl.LOCK_EX)
of2.write(sys.argv[2]+' v_fam = '+str(hits)+'\n')
fcntl.flock(of2, fcntl.LOCK_UN)
of2.close()