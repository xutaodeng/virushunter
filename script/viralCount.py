#!/usr/bin/env python
import sys
from operator import itemgetter, attrgetter
import linecache
from collections import defaultdict

def getBarcode(line, cwd):
	# line='>Contig1_170412_charlys_metagenomics_2_Dani12Torque_teno_virus_13 length=798'
	# cwd='170412_charlys_metagenomics_2'
	# virname='Torque_teno_virus_13'
	parts=line.strip().split()
	# virname=parts[1]
	readid=parts[0]
	barcode=readid[(readid.find(cwd)+1+len(cwd)):]
	return barcode

def ViralCount(infile, outfile, base, virname):
	counts=defaultdict(int)
	f = open(infile, 'r')
	of = open(outfile, 'a')
	for line in f:
		if line.strip().startswith('>'):
			barcode =getBarcode(line, base)
			counts[barcode]+=1
	f.close()
	for barcode in counts.keys():
		of.write(barcode+'\t'+virname+'\t'+str(counts[barcode])+'\n')
	of.close()

if __name__ == '__main__':
	infile, outfile, cwd, virname= sys.argv[1:5]
	ViralCount(infile, outfile, cwd, virname)

