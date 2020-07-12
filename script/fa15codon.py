#!/usr/bin/env python
import sys
f=open(sys.argv[1], 'r')
of=open(sys.argv[2],'w')

seq=[]
id=''
good, bad=0,0
for line in f:
	if line.strip().startswith('>'):
		id1=line.strip()
		#seq1='\n'.join(seq)
		seq1=''.join(seq)
		if len(seq1)>100:
			seq2=seq1[45:]
			print >>of, id
			print >>of, seq2
			if seq1[0:3]=='ATG':
				good+=1
			else:
				bad+=1
		id=id1
		seq=[]
	else:
		seq.append(line.strip())

seq1=''.join(seq)
if len(seq1)>100:
	seq2=seq1[45:]
	print >>of, id
	print >>of, seq2

print 'good', good, 'bad', bad

f.close()
of.close()
