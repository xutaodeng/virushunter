#!/usr/bin/env python
import sys, os
import os.path
from operator import itemgetter, attrgetter
from collections import defaultdict

def mergeTable(infiles, outfile, index):
	dic={}
	of=open(outfile, 'w')
	of.write('cat\tfamily\tgenus\tspecies')
	for infile in infiles:
		file=os.path.basename(infile).split('/')[-1][0:(-17)]
		of.write('\t'+file)
	of.write('\n')
	for infile in infiles:
		print 'reading...', infile
		file=os.path.basename(infile).split('/')[-1][0:(-17)]
		f=open(infile, 'r')
		for line in f:
			parts=line.strip().split()
			key='\t'.join(parts[0:4])
			if not dic.has_key(key):
				dic[key]={}
			dic[key][file]=parts[index]
		f.close()
	for key in dic.keys():
		of.write(key)
		for infile in infiles:
			file=os.path.basename(infile).split('/')[-1][0:(-17)]
			if dic[key].has_key(file):
				of.write('\t'+dic[key][file])
			else:
				of.write('\t0')
		of.write('\n')
	of.close()

if __name__ == '__main__':
	infiles, outfile = sys.argv[1:-1], sys.argv[-1]
	mergeTable(infiles, outfile+'_e-10.xls', -1)
	mergeTable(infiles, outfile+'_e-5.xls', -2)
	mergeTable(infiles, outfile+'_e-2.xls', -3)
