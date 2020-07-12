#!/usr/bin/env python
import sys
import os

def fq_pair_clean(f1, f2, of1, of2, length=1):
	i=1
	while 1:
		id=f1.readline().strip()
		if id=='': break
		seq=f1.readline().strip()
		qid=f1.readline().strip()
		qual=f1.readline().strip()
		id2=f2.readline().strip()
		seq2=f2.readline().strip()
		qid2=f2.readline().strip()
		qual2=f2.readline().strip()
		if len(seq)>length or len(seq2)>length:
			of1.write('\n'.join([id, seq, qid, qual])+'\n')
			of2.write('\n'.join([id2, seq2, qid2, qual2])+'\n')
def fq_single_clean(f1, of1, length=1):
	i=1
	while 1:
		id=f1.readline().strip()
		if id=='': break
		seq=f1.readline().strip()
		qid=f1.readline().strip()
		qual=f1.readline().strip()
		if len(seq)>length:
			of1.write('\n'.join([id, seq, qid, qual])+'\n')

if __name__ == "__main__":
	if len(sys.argv)>=5:
		f1=open(sys.argv[1], 'r')
		f2=open(sys.argv[2], 'r')
		of1=open(sys.argv[3], 'w')
		of2=open(sys.argv[4], 'w')
		length=int(sys.argv[5])
		fq_pair_clean(f1, f2, of1, of2, length)
		f1.close()
		f2.close()
		of1.close()
		of2.close()
	else:
		f1=open(sys.argv[1], 'r')
		of1=open(sys.argv[2], 'w')
		length=int(sys.argv[3])
		fq_single_clean(f1,of1, length)
		f1.close()
		of1.close()
