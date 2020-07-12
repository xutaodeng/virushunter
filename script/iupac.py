#!/usr/bin/env python
import sys
IUPAC={
'Y': 'C', 
'S': 'G',
'K': 'G', 
'B': 'C' }

nuc=set(['A','C','G','T'])

f=open(sys.argv[1], 'r')
of=open(sys.argv[2],'w')
i=0
for line in f:
	if line.strip().startswith('>') :
		id='_'.join(line.strip()[1:].split())
		print >>of, '>', id
	else:
		seq=line.strip()
		ns=[]
		for s in seq:
			if s in nuc:
				ns.append(s)
			elif IUPAC.has_key(s):
				ns.append(IUPAC[s])
			else:
				ns.append('A')
		print >>of, ''.join(ns)
f.close()
of.close()