#!/usr/bin/env python
import subprocess
import time
import sys
import os
from operator import itemgetter

#one mutation per transcript per line summarization
def readVCF(vcf, mode):
	f=open(vcf, 'r')
	VCF={}
	dps=[]
	for line in f:
		if line.startswith('#'): continue
		parts = line.strip().split()
		try: chro, pos, rs, ref, alt, qual, status, info, format, data=parts
		except: continue
		if mode == 'ALL':
			pass
			#VCF[chro+'\t'+pos]=(ref, alt)
		else: 
			iparts=info.split(';')
			#dp=iparts[6].split('=')[1]
			try:
				dp2=data.split(':')[2]
				dps.append(int(dp2))
			except: pass
			#print 'dp', dp, 'dp2', dp2
			#assert (dp==dp2)
			try: 
				if int(dp2)<10: continue
			except: pass
			if status!='PASS': continue
			eff=iparts[-1][4:]
			if 'HIGH' in eff or 'MODERATE' in eff:
				effparts=eff.split(',')
				annot=[]
				for part in effparts:
					if ('MODERATE' in part) or ('HIGH' in part):
						#annot.append(part)
						#annot=','.join(annot)
						VCF[chro+'\t'+pos+'\t'+ref+'\t'+alt+'\t'+part]=(ref, alt, part)
	print 'average depth=', sum(dps)/len(dps)
	return VCF

# DNA=["BF","CV","DT","HH","HKY","JC"]
# #RNA=["Pt_1","Pt_2","Xeno_1","Xeno_2","Xeno_3_1st","Xeno_3_2nd","Xeno_4","Xeno_5","Xeno_6","Xeno_7"]
# RNA=["Xeno_1","Xeno_2","Xeno_3_1st","Xeno_3_2nd","Xeno_4","Xeno_5","Xeno_6","Xeno_7"]
# # s=DNA+RNA

# # ss=["/mnt/san2/deng2/RNA_Seq_James/VCFs/"+x+'.vcf' for x in s]

# s=RNA
# ss=["/mnt/cluster2/deng2/RNASeq_James/"+x+'.vcf' for x in s]

ss=sys.argv[1:]

print ss
s=[x.split('/')[-3] for x in ss]
rvals=set([])
allrval={}
for sample in ss:
	rval=readVCF(sample, 'IMPORTANT')
	rvals=set(rval.keys())|rvals
	allrval[sample]=rval
of=open('vcfsummary.txt', 'w')

allout=[]
out=['chro', 'pos', 'ref', 'alt', 'Effect']
for sample in s:
	out.append(sample)
out.append('Count')
of.write('\t'.join(out)+'\n')
for key in rvals:
	out=[key]
#	found=False
	count=0
	for sample in ss:
		if allrval[sample].has_key(key):
			count+=1
			# if not found: 
				# ref, alt, eff = allrval[sample][key]
				# #out.append(ref)
				# #out.append(alt)
				# found=True
			out.append('1')
		else:
			out.append('0')
	#out.append(eff)
	allout.append((count, out))
allout.sort(key=itemgetter(0),reverse=True)
for (count, out) in allout:
	of.write('\t'.join(out)+'\t'+str(count)+'\n')
of.close()


