#!/usr/bin/env python
import sys
import os, os.path
import gzip
filename=sys.argv[1]
print filename
if filename.endswith('.gz'):f = gzip.open(filename, 'rb')
else: f=open(filename, 'r')
fileID = sys.argv[2]
#of=gzip.open(sys.argv[3],'ab')
of=open(sys.argv[3],'w')
#fileID=os.path.basename(sys.argv[2]).rsplit('.',1)[0]

i=0 #this has to be consistent with blast_trim.py
for line in f:
	i+=1
	if i%4==2:
		seq=line.strip()
		if len(seq)>=10:
			of.write('>'+fileID+'_'+str(i/4)+'\n')
			of.write(seq+'\n')
f.close()
of.close()

