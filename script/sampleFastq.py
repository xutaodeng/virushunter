#!/usr/bin/env python
import sys
import random
import gzip

f1=gzip.open(sys.argv[1], 'rb')
f2=gzip.open(sys.argv[2], 'rb')
of1=open(sys.argv[3], 'w')
of2=open(sys.argv[4],'w')
nreads=int(sys.argv[5])
i=0
for line in f1:
	i+=1
f1.close()
treads=i/4
if treads < nreads: 
	print 'treads', treads, 'nreads', nreads, sys.argv[1]
	sys.exit()
ss = set(random.sample(xrange(treads), nreads))

f1=gzip.open(sys.argv[1], 'rb')
i=0
k=0
for line in f1:
	line2 = f2.readline()
	i+=1
	if i%4==1:
		id=line.strip()
		id2=line2.strip()
	elif i%4==2:
		seq=line.strip()
		seq2=line2.strip()
	elif i%4==3:
		qid=line.strip()
		qid2=line2.strip()
	elif i%4==0:
		qseq=line.strip()
		qseq2=line2.strip()
		if i/4-1 in ss:
			print >>of1, '\n'.join([id, seq, qid, qseq])
			print >>of2, '\n'.join([id2, seq2, qid2, qseq2])
			k+=1

print 'i=', i, 'k=', k
f1.close()
f2.close()
of1.close()
of2.close()


