#!/usr/bin/env python
import os
import math
import sys

f=open(sys.argv[1], 'r') #'peak.bed', 'r')
of=open(sys.argv[2], 'w') #'enhancer.out', 'w')

#chr1    1250728 1250729 MACS_peak_13    14.00

enhChro, enhStart, enhEnd, enhName='', 0, 0, []
chro, start, end='', 0, 0
intensity=0
i=0
for line in f:
	parts = line.strip().split()
	chro, start, end, pname, height = parts
	start, end, height = int(start), int(end), float(height)
	
	if i==0: #first line initizlize
		enhChro, enhStart, enhEnd, intensity, enhName=chro, start, end, (end-start)*height,[pname]
		print 
	elif chro==enhChro and math.fabs(start-enhEnd) <12500: #connect
		enhEnd=end
		enhName.append(pname)
		intensity+=(end-start)*height
	else: #not connect, output the enhancer
		of.write(enhChro+'\t'+str(enhStart)+'\t'+str(enhEnd)+'\t'+str(intensity)+'\t'+','.join(enhName)+'\n')
		enhChro, enhStart, enhEnd, intensity, enhName=chro, start, end, (end-start)*height,[pname]
	i+=1

of.write(enhChro+'\t'+str(enhStart)+'\t'+str(enhEnd)+'\t'+str(intensity)+'\t'+','.join(enhName)+'\n')

f.close()
of.close()

