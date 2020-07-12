#!/usr/bin/env python
import sys
import gzip
import random

if sys.argv[1].endswith('.gz'):
	f = gzip.open(sys.argv[1], 'rb')
else:
	f=open(sys.argv[1], 'r')

if sys.argv[2].endswith('.gz'):
	of = gzip.open(sys.argv[2], 'wb')
else:
	of=open(sys.argv[2],'w')

nreads=int(sys.argv[3]) #should be 2million

i=0
for line in f:
	i+=1
f.close()
total_reads =i/4 #should be around 8million
random.seed()
sel=set(random.sample(xrange(total_reads), nreads))
 
if sys.argv[1].endswith('.gz'):
	f = gzip.open(sys.argv[1], 'rb')
else:
	f=open(sys.argv[1], 'r')
i=0
for line in f:
	i+=1
	if i%4==1:
		id=line.strip()
	elif i%4==2:
		seq=line.strip()
	elif i%4==3:
		qid=line.strip()
	elif i%4==0:
		qseq=line.strip()
		if i/4 in sel:
			of.write('\n'.join([id,seq,qid,qseq])+'\n')
f.close()
of.close()

