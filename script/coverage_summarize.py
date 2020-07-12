#!/usr/bin/env python
import sys
import gzip

f1 = open(sys.argv[1], 'r')
f2 = open(sys.argv[2], 'r')
f3 = open(sys.argv[3], 'r')
of=open('out.txt', 'w')
#sample /mnt/san2/HIV/160208_Leila_3run/9024-11.sam unaligned 165842 aligned reads 189 aligned percentage 0.1 %
d1,d2,d3={},{},{}
for line in f1:
	parts=line.strip().split()
	sample=parts[1].split('/')[-1].split('.')[0]
	depth=float(parts[-10])
	unaligned=int(parts[-8])
	aligned=int(parts[-5])
	per=float(parts[-2])
	d1[sample]=(unaligned, aligned, per, depth)

for line in f2:
	parts=line.strip().split()
	sample=parts[1].split('/')[-1].split('.')[0]
	depth=float(parts[-10])
	unaligned=int(parts[-8])
	aligned=int(parts[-5])
	per=float(parts[-2])
	d2[sample]=(unaligned, aligned, per, depth)
	
for line in f3:
	parts=line.strip().split()
	try: sample=parts[1].split('/')[-1].split('.')[0]
	except: print line; sys.exit()
	depth=float(parts[-10])
	unaligned=int(parts[-8])
	aligned=int(parts[-5])
	per=float(parts[-2])
	d3[sample]=(unaligned, aligned, per, depth)

for k in d1.keys():
	print k, d1[k][2], d2[k][2], d3[k][2]
	print >> of, k,  d1[k][1], d2[k][1], d3[k][1], d1[k][2],  d2[k][2], d3[k][2], d1[k][3],  d2[k][3], d3[k][3]
f1.close()
f2.close()
of.close()

