#!/usr/bin/env python
import sys
def phredfix(qseq):
	return ''.join([chr(ord(x)-33+64) for x in qseq])
f=open(sys.argv[1], 'r')
of=open(sys.argv[2],'w')

i=0
for line in f:
	i+=1
	if i%4==0:
		print >>of, phredfix(line.strip())
	else:
		print >>of, line.strip()
		
f.close()
of.close()

