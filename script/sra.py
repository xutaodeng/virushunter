#!/usr/bin/env python
import sys, gzip
infile = sys.argv[1]
if infile.endswith('.gz'):
	f=gzip.open(sys.argv[1], 'rb')
else:
	f=open(sys.argv[1], 'r')
of=gzip.open(sys.argv[3],'ab')
id=sys.argv[2]

i=0
for line in f:
	i+=1
	if line.strip()=='': continue
	if i%4==1:
		of.write(line.strip().split()[0]+' '+id+'\n')
	else:
		of.write(line.strip()+'\n')
if i%4!=0:
	print 'invalid fastq', i%4, infile
f.close()
of.close()

# print 'checking sra fastq...'
# f=gzip.open('sra.fq.gz','rb')
# i=0
# for line in f:
	# i+=1
	# if i%4==1 and not line.strip().startswith('@'): print 'wrong header', line
	# elif i%4==3 and not line.strip().startswith('+'): print 'wrong qheader', line
	# elif i%4==2: sl=len(line.strip())
	# elif i%4==0: 
		# ql=len(line.strip())
		# assert ql==sl
# if i%4!=0:
	# print 'invalid sra fastq'
# f.close()