#!/usr/bin/env python

from collections import defaultdict
import operator
import sys
import os
import re

def fixfa(infile, outfile): #illumina 33 or 64
	f = open(infile, 'r')
	of = open(outfile, 'w')
	of.write(f.readline())#header
	seq=[]
	for line in f:
		seq.append(line.strip())
	seq=''.join(seq)
	nrow = len(seq)/80
	if len(seq)%80 !=0: nrow+=1
	for i in xrange(nrow):
		start=i*80
		end=i*80+80
		try: of.write(seq[start:end]+'\n')
		except: of.write(seq[start:]+'\n')
	f.close()
	of.close()

if __name__ == '__main__':
	infile = sys.argv[1]
	outfile =sys.argv[2] 
	fixfa(infile, outfile)
