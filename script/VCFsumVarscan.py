#!/usr/bin/env python
import subprocess
import time
import sys
import os
from operator import itemgetter

#one mutation per line summarization
# each line is a unique genomic mutation location
def readVCF_unique(vcf, mode):
	f=open(vcf, 'r')
	VCF=[]
	dps=[]
	nnn=0
	for line in f:
		if line.startswith('#'): continue
		parts = line.strip().split()
		GTs,GQs,SDPs,DPs,RDs,ADs,FREQs,PVALs,RBQs,ABQs,RDFs,RDRs,ADFs,ADRs= \
			[],[],[],[],[],[],[],[],[],[],[],[],[],[]
		chro, pos, rs, ref, alt, qual, status, info, format=parts[0:9]
		gt = parts[9:]
		for g in gt:
			gi=g.split(':')
			if len(gi)>5:
				GT,GQ,SDP,DP,RD,AD,FREQ,PVAL,RBQ,ABQ,RDF,RDR,ADF,ADR=gi
			else:
				GT,GQ,SDP=gi
				DP,RD,AD,FREQ,PVAL,RBQ,ABQ,RDF,RDR,ADF,ADR = ['.']*11
			if GT=='./.': GT='.'
			elif GT=='0/0': GT='0'
			elif GT=='0/1': GT='1'
			elif GT=='1/1': GT='2'
			elif GT=='1/0': GT='1'
			
			GTs.append(GT)
			GQs.append(GQ)
			SDPs.append(SDP)
			DPs.append(DP)
			RDs.append(RD)
			ADs.append(AD)
			FREQs.append(FREQ)
			PVALs.append(PVAL)
			RBQs.append(RBQ)
			ABQs.append(ABQ)
			RDFs.append(RDF)
			RDRs.append(RDR)
			ADFs.append(ADF)
			ADRs.append(ADR)
		GTs='\t'.join(GTs)
		GQs='\t'.join(GQs)
		SDPs='\t'.join(SDPs)
		DPs='\t'.join(DPs)
		RDs='\t'.join(RDs)
		ADs='\t'.join(ADs)
		FREQs='\t'.join(FREQs)
		PVALs='\t'.join(PVALs)
		RBQs='\t'.join(RBQs)
		ABQs='\t'.join(ABQs)
		RDFs='\t'.join(RDFs)
		RDRs='\t'.join(RDRs)
		ADFs='\t'.join(ADFs)
		ADRs='\t'.join(ADRs)
		#except: continue
		if mode == 'ALL':
			pass
			#VCF[chro+'\t'+pos]=(ref, alt)
		else: 
			iparts=info.split(';')
			ADP=iparts[0].split('=')[1]
			WT=iparts[1].split('=')[1]
			HET=iparts[2].split('=')[1]
			HOM=iparts[3].split('=')[1]
			NC=iparts[4].split('=')[1]
			mutant=int(HET)+int(HOM)
			#dp=iparts[6].split('=')[1]
			# try:
				# dp2=info.split(';')[0].split('=')[1]
				# dps.append(int(dp2))
			# except: pass
			#print 'dp', dp, 'dp2', dp2
			#assert (dp==dp2)
			# try: 
				# if int(dp2)<10: continue
			# except: pass
			if status!='PASS': continue
			eff=iparts[-1][4:]
			
			if 'HIGH' in eff or 'MODERATE' in eff:
				part=eff.split(',')[0]
				# effparts=eff.split(',')
				# annot=[]
				# for part in effparts:
					# if ('MODERATE' in part) or ('HIGH' in part):
				eeff=part.split('(')[0]
				genee=part.split('|')[5]
				#annot.append(part)
				#annot=','.join(annot) \
				ddd=(chro,pos,rs,ref,alt,eff,eeff, genee,ADP,WT,HET,HOM,NC,\
					GTs,GQs,SDPs,DPs,RDs,ADs,FREQs,PVALs,RBQs,ABQs,RDFs,RDRs,ADFs,ADRs,mutant)
				nnn=len(ddd)
				VCF.append(ddd)

	VCF.sort(key=itemgetter(nnn-1),reverse=True) #sort by number of mutants
	#print 'average depth=', sum(dps)/len(dps)
	return VCF

#per transcript per mutation
#each mutation is listed multiple times for each transcript and effect
def readVCF(vcf, mode):
	f=open(vcf, 'r')
	VCF=[]
	dps=[]
	nnn=0
	for line in f:
		if line.startswith('#'): continue
		parts = line.strip().split()
		GTs,GQs,SDPs,DPs,RDs,ADs,FREQs,PVALs,RBQs,ABQs,RDFs,RDRs,ADFs,ADRs= \
			[],[],[],[],[],[],[],[],[],[],[],[],[],[]
		chro, pos, rs, ref, alt, qual, status, info, format=parts[0:9]
		gt = parts[9:]
		for g in gt:
			gi=g.split(':')
			if len(gi)>5:
				GT,GQ,SDP,DP,RD,AD,FREQ,PVAL,RBQ,ABQ,RDF,RDR,ADF,ADR=gi
			else:
				GT,GQ,SDP=gi
				DP,RD,AD,FREQ,PVAL,RBQ,ABQ,RDF,RDR,ADF,ADR = ['.']*11
			if GT=='./.': GT='.'
			elif GT=='0/0': GT='0'
			elif GT=='0/1': GT='1'
			elif GT=='1/1': GT='2'
			elif GT=='1/0': GT='1'
			
			GTs.append(GT)
			GQs.append(GQ)
			SDPs.append(SDP)
			DPs.append(DP)
			RDs.append(RD)
			ADs.append(AD)
			FREQs.append(FREQ)
			PVALs.append(PVAL)
			RBQs.append(RBQ)
			ABQs.append(ABQ)
			RDFs.append(RDF)
			RDRs.append(RDR)
			ADFs.append(ADF)
			ADRs.append(ADR)
		GTs='\t'.join(GTs)
		GQs='\t'.join(GQs)
		SDPs='\t'.join(SDPs)
		DPs='\t'.join(DPs)
		RDs='\t'.join(RDs)
		ADs='\t'.join(ADs)
		FREQs='\t'.join(FREQs)
		PVALs='\t'.join(PVALs)
		RBQs='\t'.join(RBQs)
		ABQs='\t'.join(ABQs)
		RDFs='\t'.join(RDFs)
		RDRs='\t'.join(RDRs)
		ADFs='\t'.join(ADFs)
		ADRs='\t'.join(ADRs)
		#except: continue
		if mode == 'ALL':
			pass
			#VCF[chro+'\t'+pos]=(ref, alt)
		else: 
			iparts=info.split(';')
			ADP=iparts[0].split('=')[1]
			WT=iparts[1].split('=')[1]
			HET=iparts[2].split('=')[1]
			HOM=iparts[3].split('=')[1]
			NC=iparts[4].split('=')[1]
			mutant=int(HET)+int(HOM)
			#dp=iparts[6].split('=')[1]
			# try:
				# dp2=info.split(';')[0].split('=')[1]
				# dps.append(int(dp2))
			# except: pass
			#print 'dp', dp, 'dp2', dp2
			#assert (dp==dp2)
			# try: 
				# if int(dp2)<10: continue
			# except: pass
			if status!='PASS': continue
			eff=iparts[-1][4:]
			
			if 'HIGH' in eff or 'MODERATE' in eff:
				effparts=eff.split(',')
				annot=[]
				for part in effparts:
					if ('MODERATE' in part) or ('HIGH' in part):
						eeff=part.split('(')[0]
						genee=part.split('|')[5]
						#annot.append(part)
						#annot=','.join(annot) \
						ddd=(chro,pos,rs,ref,alt,part,eeff, genee,ADP,WT,HET,HOM,NC,\
							GTs,GQs,SDPs,DPs,RDs,ADs,FREQs,PVALs,RBQs,ABQs,RDFs,RDRs,ADFs,ADRs,mutant)
						nnn=len(ddd)
						VCF.append(ddd)
	
	VCF.sort(key=itemgetter(nnn-1),reverse=True) #sort by number of mutants
	#print 'average depth=', sum(dps)/len(dps)
	return VCF
	
# DNA=["BF","CV","DT","HH","HKY","JC"]
# #RNA=["Pt_1","Pt_2","Xeno_1","Xeno_2","Xeno_3_1st","Xeno_3_2nd","Xeno_4","Xeno_5","Xeno_6","Xeno_7"]
# RNA=["Xeno_1","Xeno_2","Xeno_3_1st","Xeno_3_2nd","Xeno_4","Xeno_5","Xeno_6","Xeno_7"]
# # s=DNA+RNA

# # ss=["/mnt/san2/deng2/RNA_Seq_James/VCFs/"+x+'.vcf' for x in s]

# s=RNA
# ss=["/mnt/cluster2/deng2/RNASeq_James/"+x+'.vcf' for x in s]

def readMap():
	f=open('map.txt', 'r')
	ma={}
	for line in f:
		print line
		s, n = line.strip().split()
		ma[s]=n
	return ma

ss=sys.argv[1:] #sample names

try: ma=readMap()
except:
	ma={}
	for sample in ss:
		ma[sample]=sample

ns=[]
for sample in ss:
	ns.append(ma[sample])
ss=ns
print len(ss), 'ss', ss

rval=readVCF_unique('all_rs_snp.vcf', 'IMPORTANT')
of=open('vcfsummarySNP_mutperLine.txt', 'w')

allout=[]
out=['chro', 'pos', 'rs', 'ref', 'alt', 'Effect', 'Type', 'Gene', 'ADP', 'WT', 'HET', 'HOM', 'NC']

for sample in ss:
	out.append('GT:'+sample)
for sample in ss:
	out.append('GQ:'+sample)
for sample in ss:
	out.append('SDP:'+sample)
for sample in ss:
	out.append('DP:'+sample)
for sample in ss:
	out.append('RD:'+sample)
for sample in ss:
	out.append('AD:'+sample)
for sample in ss:
	out.append('FREQ:'+sample)
for sample in ss:
	out.append('PVAL:'+sample)
for sample in ss:
	out.append('RBQ:'+sample)
for sample in ss:
	out.append('ABQ:'+sample)
for sample in ss:
	out.append('RDF:'+sample)
for sample in ss:
	out.append('RDR:'+sample)
for sample in ss:
	out.append('ADF:'+sample)
for sample in ss:
	out.append('ADR:'+sample)
# GTs,GQs,SDPs,DPs,RDs,ADs,FREQs,PVALs,RBQs,ABQs,RDFs,RDRs,ADFs,ADRs
out.append('MutCount')
of.write('\t'.join(out)+'\n') #header
for item in rval:
	ii=[str(x) for x in item]
	of.write('\t'.join(ii)+'\n')
of.close()

#rval=readVCF('all_rs_indel.vcf', 'IMPORTANT')
rval=readVCF_unique('all_rs_indel.vcf', 'IMPORTANT')

of=open('vcfsummaryIndel_mutperLine.txt', 'w')

allout=[]
out=['chro', 'pos', 'rs', 'ref', 'alt', 'Effect', 'Type', 'Gene', 'ADP', 'WT', 'HET', 'HOM', 'NC']

for sample in ss:
	out.append('GT:'+sample)
for sample in ss:
	out.append('GQ:'+sample)
for sample in ss:
	out.append('SDP:'+sample)
for sample in ss:
	out.append('DP:'+sample)
for sample in ss:
	out.append('RD:'+sample)
for sample in ss:
	out.append('AD:'+sample)
for sample in ss:
	out.append('FREQ:'+sample)
for sample in ss:
	out.append('PVAL:'+sample)
for sample in ss:
	out.append('RBQ:'+sample)
for sample in ss:
	out.append('ABQ:'+sample)
for sample in ss:
	out.append('RDF:'+sample)
for sample in ss:
	out.append('RDR:'+sample)
for sample in ss:
	out.append('ADF:'+sample)
for sample in ss:
	out.append('ADR:'+sample)
# GTs,GQs,SDPs,DPs,RDs,ADs,FREQs,PVALs,RBQs,ABQs,RDFs,RDRs,ADFs,ADRs
out.append('MutCount')
of.write('\t'.join(out)+'\n') #header
for item in rval:
	ii=[str(x) for x in item]
	of.write('\t'.join(ii)+'\n')
of.close()

