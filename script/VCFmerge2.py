#!/usr/bin/env python
import subprocess
import time
import sys
import os
from operator import itemgetter

#one mutation per transcript per line summarization
def readVCFSummary(sumfile):
	f=open(sumfile, 'r')
	header=f.readline().strip().split()
	ns = (len(header)-14)/14 #number of samples
	header='\t'.join(header[13:(13+ns)])
	rval={}
	for line in f:
		parts = line.strip().split()
		chro, pos,rs,ref, alt, Effect, Type, Gene, ADP, WT, HET, HOM, NC =parts[0:13]
		geno='\t'.join(parts[13:(13+ns)])
		mutCount= int(parts[-1])
		WT, HET, HOM, NC = int(WT), int(HET), int(HOM), int(NC)
		key=chro+'_'+pos+'_'+Effect
		rval[key]=(geno, mutCount, chro, pos,rs,ref, alt, Effect, Type, Gene, ADP, WT, HET, HOM, NC)
	print 'num mut', len(rval)
	print 'num samp', ns
	return (ns, header, rval)

def Merge(sum1, sum2, out):
	ns1, header1, r1= sum1
	ns2, header2, r2=sum2
	of=open(out, 'w')
	key1=r1.keys()
	key2=r2.keys()
	ck=set(key1) | set(key2)
	ak=set(key1) & set(key2)
	rval=[]
	header= ['mutCount', 'chro', 'pos','rs','ref', 'alt', 'Effect', 'Type', 'Gene', 'WT', 'HET', 'HOM', 'NC', header1, header2]
	of.write('\t'.join(header)+'\n')
	for key in ck:
		if r1.has_key(key):
			geno1, mutCount1, chro, pos,rs,ref, alt, Effect, Type, Gene, ADP1, WT1, HET1, HOM1, NC1 = r1[key]
		else:
			geno1, mutCount1, WT1, HET1, HOM1, NC1 = '\t'.join(ns1*['*']), 0, 0,0,0,0

		if r2.has_key(key):
			geno2, mutCount2, chro, pos,rs,ref, alt, Effect, Type, Gene, ADP2, WT2, HET2, HOM2, NC2 = r2[key]
		else:
			geno2, mutCount2, WT2, HET2, HOM2, NC2 = '\t'.join(ns2*['*']), 0, 0,0,0,0
		rval.append([(mutCount1+mutCount2)*1000000000000000+int(pos), str(mutCount1+mutCount2), chro, pos,rs,ref, alt, Effect, Type, Gene, str(WT1+WT2), str(HET1+HET2), str(HOM1+HOM2), str(NC1+NC2), geno1, geno2])
	print 'total variants', len(ck)
	print 'common variants', len(ak)
	rval.sort(key=itemgetter(0),reverse=True) #sort by number of mutants
	#print 'average depth=', sum(dps)/len(dps)
	for item in rval:
		ii=item[1:]
		of.write('\t'.join(ii)+'\n')
	of.close()


sum1 = readVCFSummary(sys.argv[1]) #5
sum2 = readVCFSummary(sys.argv[2]) #6
Merge(sum1, sum2, sys.argv[3])
