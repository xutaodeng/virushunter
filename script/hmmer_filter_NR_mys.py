#!/usr/bin/env python
from collections import defaultdict
import operator
import sys
import os
import os.path
import linecache
import re

def FilterLines(submyspre, submys, hits): 
	f = open(submyspre, 'r')
	of = open(submys, 'w')
	nmys=0
	found=0
	out=False
	for line in f:
		if line.strip().startswith('>'):
			header = line.strip()[1:].strip()
			if header not in hits:
				of.write(line.strip()+'\n')
				nmys+=1
				out=True
			else: found+=1
		elif out:
			of.write(line.strip()+'\n')
			out=False
	print submys, 'n_mys =', nmys, 'n_found =', found


def hmmer_parser(fname):
	f=open(fname, 'r')
	start=False
	queryset=set()
	for line in f:
		if line.strip().startswith ("Query:"):
			queryprot=line.strip().split()[1]
			querydna = queryprot.rsplit('_', 1)[0][:-4]
			start=False
		if not start and line.strip().startswith(">>"):
			start=True
			queryset.add(querydna)
	f.close()
	return queryset

if __name__ == '__main__': 
	submysxmlpre=sys.argv[1] #xml pre mys
	submyspre = sys.argv[2] #pre mys fasta
	submys = sys.argv[3] #mys fasta filtered
	hits = hmmer_parser(submysxmlpre)
	FilterLines(submyspre, submys, hits)