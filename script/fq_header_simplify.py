#!/usr/bin/env python
#clean header of fa files, only keep the first element
import sys
f=open(sys.argv[1], 'r')
of=open(sys.argv[2], 'w')
i=0
for line in f:
	i+=1
	if line.strip().startswith('>'):
		print >>of, '>s'+str(i) 
	else:
		print >>of, line.strip()
f.close()
of.close()
