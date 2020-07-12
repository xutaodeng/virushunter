#!/usr/bin/env python
import sys
f=open(sys.argv[1], 'r')
of=open(sys.argv[2],'w')
i=0
seq=[]
for line in f:
	if line.strip().startswith('>'):
		i+=1
		id='@seq'+str(i)
		seq1=''.join(seq)
		if seq1!='':
			print >>of, id
			print >>of, seq1
			print >>of, '+'
			print >>of, ''.join(len(seq1)*['I'])
			seq=[]
	else:
		seq.append(line.strip())

print >>of, id
print >>of, seq1
print >>of, '+'
print >>of, ''.join(len(seq1)*['I'])

f.close()
of.close()