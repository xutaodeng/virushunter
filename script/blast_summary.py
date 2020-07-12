#!/usr/bin/env python

from collections import defaultdict
import operator
import sys
import os
import os.path
import linecache
import re, string

def CacheLines(input, output): 
	VE={}#defaultdict(int)
	vs=set([])
	f = open(input, 'r')
	of=open(output, 'w')
	i=0
	for line in f:
		i+=1
		if i%11 == 2:
			sample = line.strip().split('_')[-1]
			if not VE.has_key(sample): VE[sample]=defaultdict(int)
		if i%11 == 3:
			query = line.strip().split()[1]
		if i%11 == 4:
			taxonomy=line.strip().split()[-1]
			try: virus=taxonomy.split(':')[0].split('$')[1]
			except: print line; sys.exit()
			vs.add(virus)
		elif i%11==6: #veval
			ve=float(line.split()[-1])
			if ve <= 0.00001:
				VE[sample][virus]+=1
	f.close()
	
	out=['virus']
	for sample in VE.keys():
		out.append(sample)
	of.write('\t'.join(out)+'\n')
	for v in vs:
		out=[v]
		for sample in VE.keys():
			out.append(str(VE[sample][v]))
		of.write('\t'.join(out)+'\n')
	f.close()
	of.close()

if __name__ == '__main__': 
	input=sys.argv[1] #input of NR filtered blast output
	output=sys.argv[2]
	CacheLines(input, output)
	#print VE

