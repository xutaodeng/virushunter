#!/usr/bin/env python
from collections import defaultdict
import operator
import sys
import os.path
import os

i1=open(sys.argv[1], 'r') #s_G1_L001_I1_001.fastq
i2=open(sys.argv[2], 'r') #s_G1_L001_I1_002.fastq
f1=open(sys.argv[3], 'r') #myseq1.fq
f2=open(sys.argv[4], 'r') #myseq2.fq
of1=open(sys.argv[5], 'w') #m1.fq
of2=open(sys.argv[6], 'w') #m2.fq

###      +HWI-ST281_0208:1:1101:1235:2108#GTAGAG/1

ids={}
i=0
for line in i1:
	i+=1
	if i%4==1:
		id, pair = line.strip().split()
	if i%4==2:
		barcode = line.strip()
		ids[id]=barcode
i1.close()
i=0
for line in i2:
	i+=1
	if i%4==1:
		id, pair = line.strip().split()
	if i%4==2:
		barcode = line.strip()
		ids[id]=barcode
i2.close()

i=0
for line in f1:
	i+=1
	if i%4==1:
		id, pair = line.strip().split()
		pair=pair.split(':')[0]
		barcode= ids[id]
		newline=id+'#'+barcode+'/'+pair
		print >>of1, newline.strip()
	else:
		print >>of1, line.strip()
f1.close()
of1.close()

i=0
for line in f2:
	i+=1
	if i%4==1:
		id, pair = line.strip().split()
		pair=pair.split(':')[0]
		barcode= ids[id]
		newline=id+'#'+barcode+'/'+pair
		print >>of2, newline.strip()
	else:
		print >>of2, line.strip()
f2.close()
of2.close()
		