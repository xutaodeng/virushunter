#!/usr/bin/env python
import subprocess
import time
import sys
import os
from operator import itemgetter


def readSample(file):
	f=open(file, 'r')
	rval={}
	for line in f:
		parts=line.strip().split()
		chro, pos, ref, alt, eff = parts
		#rval[chro+'\t'+pos]=(ref, alt, eff)
		effparts=eff.split(',')
		annot=[]
		for part in effparts:
			if ('MODERATE' in part) or ('HIGH' in part):
				annot.append(part)
		annot=','.join(annot)
		rval[chro+'\t'+pos+'\t'+ref+'\t'+alt]=(ref, alt, annot)
	f.close()
	return rval

# s=["Xeno_6_v_BF", "Xeno_3_1st_v_CV", "Xeno_4_v_DT", "Xeno_1_v_HH", "Xeno_2_v_HKY", "Xeno_5_v_JC", "Xeno_1_v_Pt_1", "Xeno_2_v_Pt_2", "Pt_1_v_Xeno_1", "Pt_2_v_Xeno_2"]

# ss=["/mnt/san2/deng2/RNA_Seq_James/VCFs/"+x for x in s]


ss=sys.argv[1:]
s=[x.split('/')[-3] for x in ss]
print ss
rvals=set([])
allrval={}
for sample in ss:
	rval=readSample(sample)
	rvals=set(rval.keys())|rvals
	allrval[sample]=rval
of=open('diffsummary.txt', 'w')

allout=[]
out=['chro', 'pos', 'ref', 'alt']
for sample in s:
	out.append(sample)
out.append('Effect')
out.append('Count')
of.write('\t'.join(out)+'\n')
for key in rvals:
	out=[key]
	found=False
	count=0
	for sample in ss:
		if allrval[sample].has_key(key):
			count+=1
			if not found: 
				ref, alt, eff = allrval[sample][key]
				#out.append(ref)
				#out.append(alt)
				found=True
			out.append('1')
		else:
			out.append('0')
	out.append(eff)
	allout.append((count, out))
allout.sort(key=itemgetter(0),reverse=True)
for (count, out) in allout:
	of.write('\t'.join(out)+'\t'+str(count)+'\n')
of.close()

