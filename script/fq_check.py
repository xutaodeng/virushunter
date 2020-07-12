#!/usr/bin/env python
import sys

ss =set(['A','C','G','T','N'])
def testSeq(seq):
	a = set(list(seq))
	for s in a:
		assert (s in ss)

f=open(sys.argv[1], 'r')
print 'checking..', sys.argv[1]
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
		testSeq(seq)
		try: assert len(qseq) == len(seq)
		except: print len(seq), len(qseq), seq, qseq
f.close()
assert i%4==0
print sys.argv[1], 'done'
