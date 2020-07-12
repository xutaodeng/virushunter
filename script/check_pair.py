#!/usr/bin/env python
import sys
import difflib
def revcomp(seq):
	rval=[]
	for s in seq:
		if s=='A':
			rval.append('T')
		elif s=='T':
			rval.append('A')
		elif s=='G':
			rval.append('C')
		elif s=='C':
			rval.append('G')
		else:
			rval.append(s)
	rval.reverse()
	return ''.join(rval)

f1=open(sys.argv[1], 'r')
f2=open(sys.argv[2], 'r')
i=0
totalOverlap=0
totalLength=0
nr = 0 #number of overlap
for line1 in f1:
	line2 = f2.readline()
	i+=1
	if i%4==2:
		seq1=line1.strip()
		seq2=line2.strip()
		if len(seq1)<=1 or len(seq2)<=1:
			continue
		seq2=revcomp(seq2)
		m =difflib.SequenceMatcher(None, seq1, seq2).find_longest_match(0,len(seq1),0,len(seq2))
		a, b, r = m.a, m.b, m.size #a is seq1 start, b seq2 start, r is the match size
		if r<=10 or b>10 : continue
		totalOverlap+=r
		totalLength+=max(len(seq1), len(seq2))
		nr+=1
print sys.argv[1], 'ave_bp_Overlap =', totalOverlap/nr
print sys.argv[1], 'nreads_Overlap =', nr
f1.close()
f2.close()