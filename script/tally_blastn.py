#!/usr/bin/env python

from collections import defaultdict
import operator
import sys
import os
import re

def readFa(ref):
	f=open(ref, 'r')
	map={}
	for line in f:
		if line.strip().startswith('>'):
			gi=line.strip()[1:].split()[0]
			map[gi]=line.strip()
	return map

def tally(map, infiles, outfile): 
	of =open(outfile+'_2', 'w')
	of5 =open(outfile+'_5', 'w')
	of10 =open(outfile+'_10', 'w')
	counts={}
	counts5={}
	counts10={}
	for infile in infiles:
		f = open(infile, 'r')
		for line in f:
			try: 
				parts = line.strip().split()
				contig=parts[1]
				Evalue=float(parts[-2])
			except: print infile; sys.exit()
			if not counts.has_key(contig):
				counts[contig]=defaultdict(int)
				counts5[contig]=defaultdict(int)
				counts10[contig]=defaultdict(int)
			if Evalue<=1E-2: counts[contig][infile]+=1
			if Evalue<=1E-5: counts5[contig][infile]+=1
			if Evalue<=1E-10: counts10[contig][infile]+=1
			
		f.close()
	out, out5, out10=['ref'], ['ref'], ['ref']
	for infile in infiles:
		out.append(infile.split('/')[-1].split('.')[0])
		out5.append(infile.split('/')[-1].split('.')[0])
		out10.append(infile.split('/')[-1].split('.')[0])
	of.write('\t'.join(out)+'\n')
	of5.write('\t'.join(out5)+'\n')
	of10.write('\t'.join(out10)+'\n')
	for contig in counts.keys():
		out=[map[contig]]
		out5=[map[contig]]
		out10=[map[contig]]
		for infile in infiles:
			try: out.append(str(counts[contig][infile]))
			except: out.append('0')
			try: out5.append(str(counts5[contig][infile]))
			except: out5.append('0')
			try: out10.append(str(counts10[contig][infile]))
			except: out10.append('0')
			
		of.write('\t'.join(out)+'\n')
		of5.write('\t'.join(out5)+'\n')
		of10.write('\t'.join(out10)+'\n')
	of.close()
	of5.close()
	of10.close()

if __name__ == '__main__':
	ref=sys.argv[1]
	infiles = sys.argv[2:-1]
	outfile =sys.argv[-1] 
	map=readFa(ref)
	tally(map, infiles, outfile)
