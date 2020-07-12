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
key=sys.argv[2]

def readSam(sam, key):
	os.system('samtools view -c '+sam+' > 1')
	os.system('samtools view -c -F 260 '+sam+ ' > 2 ')
	of=open('read.count', 'a')

	f=open('1', 'r')
	a=f.readline().strip()
	f.close()
	f=open('2', 'r')
	b=f.readline().strip()
	f.close()
	of.write(key+' totalRead '+a+' mapped '+b+' percentage '+str(float(b)/float(a)*100)+'\n')
	of.close()

readSam(sam, key)

