#!/usr/bin/env python

import os
from collections import defaultdict
import os.path
import sys
import math

# tRNA
# transfer RNA
# Mt-tRNA
# transfer RNA located in the mitochondrial genome
# rRNA
# ribosomal RNA
# scRNA
# small cytoplasmic RNA
# snRNA
# small nuclear RNA
# snoRNA
# small nucleolar RNA
# miRNA
# microRNA precursors
# misc_RNA
# miscellaneous other RNA
# lincRNA
# Long intergenic non-coding RNAs
##################################################################
#### construct non-coding RNA database
#1. download from ensemble
#http://uswest.ensembl.org/info/data/ftp/index.html
#2. Abundant sequences from tophat download
#/mnt/cluster/mm10/Mus_musculus/UCSC/mm10/Sequence/AbundantSequences
#/mnt/cluster/hg19/Sequence/AbundantSequences
#3. cat these seuqnces cat *.fa > nc.fasta and 
#4. bowtie2-build nc.fasta nc
##################################################################
sam=sys.argv[1]
fa=sys.argv[2]
samC=sys.argv[3]

def readHeader(fa):
	dic={}
	f=open(fa, 'r')
	for line in f:
		if not line.startswith('>'): continue
		id=line.strip()[1:]
		parts=id.split()
		if len(parts)==1:
			dic[id]=id
			#>polyA
		elif id.startswith('gi|555853|gb|U13369.1|HSU13369'): #human
			dic['gi|555853|gb|U13369.1|HSU13369']='rRNA'
		elif id.startswith('gi|38176281|tpg|BK000964.1'): #mouse
			dic['gi|38176281|tpg|BK000964.1']='rRNA'
		elif id.startswith('gi'):
			dic[parts[0]] = '_'.join(parts[1:])
			#>gi|9626372|ref|NC_001422.1| Enterobacteria phage phiX174, complete genome
		else:
			#>ENST00000410344.1 ncrna:known chromosome:GRCh38:14:96384624:96384815:1 gene:ENSG00000222276.1 gene_biotype:snRNA transc
			dic[parts[0]] = parts[4].split(':')[1]
	f.close()
	print 'dic size', len(dic)
	return dic
def readSam(sam, dic, samC):
	res=defaultdict(int)
	f=open(sam, 'r')
	cc=0
	for line in f:
		cc+=1
		if cc%10000000==0: 
			sys.stdout.write('\r')
			sys.stdout.write(str(cc))
			sys.stdout.flush()
		parts=line.strip().split()
		id = parts[2]
		if id=='*': continue
		res[dic[id]]+=1
	f.close()
	of=open(samC, 'w')
	for key, val in res.items():
		of.write(key+' '+str(val)+' percentage '+str(float(val)/cc*100)+'\n')

	of.close()

dic=readHeader(fa)
readSam(sam, dic, samC)