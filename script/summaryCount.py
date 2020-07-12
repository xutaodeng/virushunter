#!/usr/bin/env python

from collections import defaultdict
import operator
import sys
import os
import os.path
import linecache
import fcntl

def process(countfiles, outfile):
	of=open(outfile, 'w')
	a=defaultdict(int)
	csv=False
	for cf in countfiles:
		#print 'cf', cf
		f=open(cf, 'r')
		for line in f:
			#print line
			if line.startswith('category,'): csv=True; continue
			if csv:
				key, count = line.strip().rsplit(',', 1)
			else:
				try:key, count = line.strip().split()
				except: print 'line', line; sys.exit()
			count= int(count)
			a[key]+=count
		f.close()
	for key in a.keys():
		if csv: of.write(key+','+str(a[key])+'\n')
		else: of.write(key+'\t'+str(a[key])+'\n')
	of.close()

if __name__ == '__main__': 
	countfiles=sys.argv[1:-1]
	outfile=sys.argv[-1]
	process(countfiles, outfile)