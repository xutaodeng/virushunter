#!/usr/bin/env python

import sys
import os
import os.path
from optparse import OptionParser
from collections import defaultdict, deque
from bisect import bisect_left
import re
from operator import itemgetter

def processPairedSAM(filename, filename2): # first scan to get the mutation positions
	un, p1, p2, paired=0,0,0,0
	if filename[-3:]==".gz":
		import gzip
		f=gzip.open(filename)
	else: f=open(filename, 'r')
	# of=open(outfile, 'w')
	if filename2[-3:]==".gz":
		import gzip
		f2=gzip.open(filename2)
	else: f2=open(filename2, 'r')
	# of2=open(outfile2, 'w')
	count=0
	for line in f:
		line2=f2.readline()
		if line[0]=='@': continue
		count+=1
		parts=line.strip().split('\t')
		parts2=line2.strip().split('\t')
		(name, flag, chro, start, mapq, cigar)=parts[0:6]
		(name2, flag2, chro2, start2, mapq2, cigar2)=parts2[0:6]
		name='@'+name
		name2='@'+name2
		#print 'cigar', cigar
		try: seq = parts[9]
		except: print line; seq='A'
		try: seq2 = parts2[9]
		except: print line2; seq2='A'
		qual=parts[10]
		qual2=parts2[10]

		if chro != '*'  and chro2!='*':
			paired+=1
		elif chro !='*':
			p1+=1
		elif chro2!='*':
			p2+=1
		else:
			un+=1
	f.close()
	print filename, 'paired', paired, 'p1Only', p1, 'p2Only', p2, 'unaligned', un 

if __name__ == "__main__":
	processPairedSAM(sys.argv[1], sys.argv[2]) #pair end

