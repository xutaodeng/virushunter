#!/usr/bin/env python
import subprocess
import time
import sys
import os
from operator import itemgetter

def readdbSNP(dbsnpfile):
	rval={}
	f=open(dbsnpfile, 'r')
	for line in f:
		if line.startswith('#'): continue
		chro, pos, rs = line.strip().split()[0:3]
		rval[chro+':'+pos]=rs
	return rval

#one mutation per transcript per line summarization
def annoVCF(vcf, vcfout, dbSNP):
	f=open(vcf, 'r')
	of=open(vcfout, 'w')
	for line in f:
		if line.startswith('#'): of.write(line); continue
		parts = line.strip().split()
		chro, pos, rs, ref, alt, qual, status, info, format=parts[0:9]
		if dbSNP.has_key(chro+':'+pos):
			parts[2] = dbSNP[chro+':'+pos]
			of.write('\t'.join(parts)+'\n')
		else:
			of.write(line)
	f.close()
	of.close()

dbsnpfile=sys.argv[1] #sample names
vcf=sys.argv[2] #sample names
vcfout=sys.argv[3] #sample names
dbSNP=readdbSNP(dbsnpfile)
annoVCF(vcf, vcfout, dbSNP)
