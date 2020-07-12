#!/usr/bin/env python
#combine two fq files but filter dups, generated from blastx and blastn
import sys
f1=open(sys.argv[1], 'r')
f2=open(sys.argv[2],'r')
i=0
ids=set()
for line in f1:
	i+=1
	if i%4==1:
		ids.add(line.strip())
	print line.strip()

i=0
for line in f2:
	i+=1
	if i%4==1:
		id=line.strip()
		if id in ids:
			dup=True
		else:
			dup=False
			print id
	elif not dup:
		print line.strip()
			
f1.close()
f2.close()