#!/usr/bin/env python

import sys
import os
import os.path
from optparse import OptionParser
from collections import defaultdict, deque
from bisect import bisect_left
import re
from operator import itemgetter


def CIGARMatch(cigar):
#first scan to get the mutation dictionary
	N=re.findall(r'[0-9]+', cigar) #numbers
	L=re.findall(r'[DNHPMIS]', cigar) #letters
	for i in range(len(L)):
		l,n = L[i], int(N[i])
		if l=='M' and int(n)>=50:
			return True
	return False

def processTrimSAM(filename, outfile): # first scan to get the mutation positions
	filter=0
	if filename[-3:]==".gz":
		import gzip
		f=gzip.open(filename)
	else: f=open(filename, 'r')
	of=open(outfile, 'w')
	for line in f:
		if line[0]=='@': continue	   
		parts=line.strip().split('\t')
		(name, flag, chro, start, mapq, cigar)=parts[0:6]
		name='@'+name
		try: seq = parts[9]
		except: seq='A'
		qual=parts[10]
		start=int(start)
		if len(seq)==1:
			filter+=1
			print >>of, '\n'.join([name, 'A', '+','A'])
		elif chro == '*':
			print >>of, '\n'.join([name, seq, '+',qual])
		elif chro in virus:
			print >>of, '\n'.join([name, seq, '+',qual])
		else:
			if CIGARMatch(cigar): #match non virus seq in NT
				print >>of, '\n'.join([name, 'A', '+','A'])
				pie[chro]+=1
			else: #match non virus but not long enough
				print >>of, '\n'.join([name, seq, '+',qual])
				pie['*']+=1
	f.close()
	print filename, 'filtered', filter
	sum=0
	for key in pie.keys():
		print key, pie[key]
		sum+=pie[key]
	print 'candidate virus percent', 100.0*(pie['unclass']+ pie['Viruses'])/sum+'%'
	

#only consider the quality of mutant bases, i.e., not '.,'
if __name__ == "__main__":
	samfile, outfile =sys.argv[1], sys.argv[2]
	processTrimSAM(samfile, outfile)