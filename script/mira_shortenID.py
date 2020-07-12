#!/usr/bin/env python
import sys
if sys.argv[1].endswith('.gz'): f = gzip.open(sys.argv[1], 'rb')
else: f=open(sys.argv[1], 'r')
of=open(sys.argv[2],'w')
pair = sys.argv[3]
i=0
for line in f:
	i+=1
	if i%4==1:
		id='@mira'+str(i/4)+' '+pair
	elif i%4==2:
		seq=line.strip()
	elif i%4==3:
		qid=line.strip()
	elif i%4==0:
		qseq=line.strip()
		of.write('\n'.join([id,seq, qid, qseq])+'\n')
f.close()
of.close()

