#!/usr/bin/env python
import os
import math
import sys
from bisect import bisect_left
from collections import defaultdict



#chr1    1612739 1615288 726905.97       MACS_peak_43,MACS_peak_44,MACS_peak_45

def readBed(f):
	print 'reading bed...'
	beds=defaultdict(list)
	starts=defaultdict(list)
	for line in f:
		parts=line.strip().split()
		chro, start, end, name, value  = parts[0:5]
		beds[chro].append([int(start), int(end), name])
		starts[chro].append(int(start))
	print 'reading bed done!'
	print beds.keys()
	return beds, starts


def readSam(f1):
	counts=defaultdict(int)
	kkk=0
	for line in f1:
		if line.startswith('@'): continue
		kkk+=1
		if kkk%1000000==0: print 'processing read ', kkk
		parts=line.strip().split()
		chro=parts[2]
		if chro=='*': continue
		if len(beds[chro])==0: continue
		bstart=int(parts[3])
		bend=bstart+30
		i=bisect_left(starts[chro], bstart)
		i=i-1
		if i<0: i=0
		s, e, name = beds[chro][i]
		if max(bstart,s )<= min(bend, e) :
			counts[chro+':'+str(s)]+=1
		#next bed
		i+=1
		if i== len(beds[chro]): continue
		
		s, e, name = beds[chro][i]
		if max(bstart,s )<= min(bend, e) :
			counts[chro+':'+str(s)]+=1
	return float(kkk)/1000000, counts

def PrintCounts2(counts1, counts2, n1, n2, of):
	for ch in beds.keys():
		for (s, e, name) in beds[ch]:
			c1=counts1[ch+':'+str(s)]
			c2=counts2[ch+':'+str(s)]
			c=c1-c2
			of.write(ch+'\t'+str(s)+'\t'+str(e)+'\t'+str(e-s)+'\t'+str(c1)+'\t'+str(c2)+'\t'+str((c1/n1-c2/n2)/(e-s)*51)+'\t'+name+'\n')

def PrintCounts1(counts1, n1, of):
	for ch in beds.keys():
		for (s, e, name) in beds[ch]:
			c1=counts1[ch+':'+str(s)]
			of.write(ch+'\t'+str(s)+'\t'+str(e)+'\t'+str(e-s)+'\t'+str(c1)+'\t'+'-'+'\t'+str((c1/n1)/(e-s)*51)+'\t'+name+'\n')

f=open(sys.argv[1], 'r') #'peak.bed', 'r')
f1=open(sys.argv[2], 'r') #'sam1', 'r', treat sam)
of=open(sys.argv[3], 'w') #'enhancer.out', 'w')

beds, starts=readBed(f)
n1, counts1=readSam(f1)

try: 
	f2=open(sys.argv[4], 'r') #'sam2', control background)
	n2, counts2=readSam(f2)
	PrintCounts2(counts1, counts2, n1, n2, of)
	f2.close()
except: 
	PrintCounts1(counts1, n1, of)

f.close()
f1.close()
of.close()
