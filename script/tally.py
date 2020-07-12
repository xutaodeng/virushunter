#!/usr/bin/env python

from collections import defaultdict
import operator
import sys
import os
import re
from operator import itemgetter, attrgetter

def tally(infiles, outfile, fas): 
	of =open(outfile, 'w')
	counts={}
	labels=[]
	try:
		for infile in infiles:
			labels.append(infile.split('/')[-1].split('.')[0])
		#try:
		lab=[(x[-1],int(x[0:-1]), y) for (x,y) in zip(labels, infiles)]
		#print lab
		lab=sorted(lab, key=itemgetter(0, 1))
		#print lab
		infiles=[x[2] for x in lab]
	except:
		pass
	
	#print infiles
	#except: pass
	for infile in infiles:
		f = open(infile, 'r')
		for line in f:
			try: contig, count = line.strip().split()
			except: print infile; print 'error'; print line; sys.exit()
			if counts.has_key(contig):
				counts[contig][infile]=count
			else:
				counts[contig]={infile:count}
		f.close()
	out=['contig']
	for infile in infiles:
		out.append(infile.split('/')[-1].split('.')[0])
	of.write('\t'.join(out)+'\n')
	print 'common', len(set(counts.keys()).intersection(set(fas)))
	for contig in fas: #counts.keys():
		#print contig
		out=[contig]
		for infile in infiles:
			try: out.append(str(counts[contig][infile]))
			except: out.append('0')
		of.write('\t'.join(out)+'\n')
	of.close()

def readFA(fa):
	fas=set([])
	f=open(fa, 'r')
	for line in f:
		if line.strip().startswith('>'):
			fas.add(line.strip().split()[0][1:])
			#chr=line.strip().rsplit('_', 1)[0][1:]
			#fas.add(chr)
	#print fas
	fas=list(fas)
	return fas

if __name__ == '__main__':
	fa=sys.argv[1]
	fas=readFA(fa)
	infiles = sys.argv[2:-1]
	outfile =sys.argv[-1] 
	tally(infiles, outfile, fas)
