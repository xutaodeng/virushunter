#!/usr/bin/env python
import sys
import gzip

if sys.argv[1].endswith('.gz'):
	f = gzip.open(sys.argv[1], 'rb')
else:
	f=open(sys.argv[1], 'r')
of=open(sys.argv[2],'w')
try:length=int(sys.argv[3])
except: length=0

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
		if len(seq) >= length:
			print >>of, '>'+id.strip()
			#print >>of, '>',id,qid,qseq
			print >>of, seq
			#print 'length', length
			#print seq
f.close()
of.close()

