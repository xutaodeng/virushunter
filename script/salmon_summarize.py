#!/usr/bin/env python

from collections import defaultdict
import operator
import sys
import os
import re

def loadAnno():
	annot={}
	f=open('/mnt/cluster/tools/Homo_sapiens.GRCh38.nc.cds.all.fa', 'r')
	for line in f:
		if line.strip().startswith('>'):
			parts=line.strip().split()
			enst=parts[0][1:]
			type=parts[1]
			info=parts[2:]
			type, symbol='',''
			for pa in info:
				if pa.startswith('gene_biotype'):
					type=pa.split(':')[1]
				if pa.startswith('gene_symbol'):
					symbol=pa.split(':')[1]
			annot[enst]=(type, symbol)
	print 'num enst', len(annot)
	return annot

def readSample():
	f=open('paste.sh', 'r')
	samples=f.readline().strip().split()[1:-2]
	print samples
	samples=[s.split('/')[0][:-4] for s in samples]
	return samples

def TransToGene(infile, outfile):
	f=open(infile, 'r')
	of=open(outfile, 'w')
	header=f.readline().strip().split()
	type, gene, enst, length, efLength = header[0:5]
	samples=header[5:]
	rval={}
	for line in f:
		parts=line.strip().split()
		type, gene, enst, length, efLength = parts[0:5]
		value=parts[5:]
		value=[float(x) for x in value]
		if not rval.has_key(gene):
			rval[gene]=len(value)*[0.0]
		for i in xrange(len(value)):
			rval[gene][i]+=value[i]
	of.write('gene\t'+'\t'.join(samples)+'\n')
	for gene in rval.keys():
		dat=rval[gene]
		dat=[str(x) for x in dat]
		of.write(gene+'\t'+'\t'.join(dat)+'\n')
	of.close()
	f.close()

def readSalmon(expfile,  annot, samples, outfile, outfile2, summaryfile):
	f=open(expfile, 'r')
	of=open(outfile, 'w')
	of2=open(outfile2, 'w')
	of3=open(summaryfile,'w')
	of.write('\t'.join(['type', 'gene', 'enst', 'length', 'efLength'])+'\t')
	of.write('\t'.join(samples)+'\n')
	of2.write('\t'.join(['type', 'gene', 'enst', 'length', 'efLength'])+'\t')
	of2.write('\t'.join(samples)+'\n')
	f.readline() #header
	# ncReads={}
	# codingReads=defaultdict(float)
	reads={}
	alltype=set()
	for line in f:
		parts=line.strip().split()
		ns=len(parts)/5
		assert ns ==len(samples)
		for i in range(ns):
			enst, length, efLength, TPM, numReads= parts[(i*5+0):(i*5+5)]
			type, gene=annot[enst]
			alltype.add(type)
			if not reads.has_key(samples[i]):
				reads[samples[i]]=defaultdict(float)
			reads[samples[i]][type]+=float(numReads)
			#codingReads[samples[i]]+=float(numReads)			
			if type =='protein_coding': 
				if i==0:
					of.write('\t'.join([type, gene, enst, length, efLength, TPM])+'\t')
					of2.write('\t'.join([type, gene, enst, length, efLength, numReads])+'\t')
				else:
					of.write(TPM+'\t')
					of2.write(numReads+'\t')
				if i==ns-1:
					of.write('\n')
					of2.write('\n')
	of3.write('sample'+'\t'+'\t'.join(alltype)+'\t'+'totalReads'+'\n')
	for samp, counts in reads.items():
		total=0
		of3.write(samp+'\t')
		for type in alltype:
			total+=reads[samp][type]
		for type in alltype:
			of3.write(str(reads[samp][type]/total*100)+'%'+'\t')
		of3.write(str(total)+'\n')

	# for key, value in ncReads.items():
		# print key, value/totalnc*100,'%'

	#print expfile, 'nc=',totalnc, 'coding=',codingReads, 'nc%=',totalnc/(codingReads+totalnc)*100,'%'

# def CompreCodon

if __name__ == '__main__': 
	annot=loadAnno()
	samples=readSample()
	readSalmon('codon.sf', annot, samples, 'codon.TPM', 'codon.reads', 'codon.summary')
	TransToGene('codon.TPM', 'codon_gene.TPM')
	TransToGene('codon.reads', 'codon_gene.reads')
	
	readSalmon('full.sf', annot, samples, 'full.TPM', 'full.reads', 'full.summary')
	TransToGene('full.TPM', 'full_gene.TPM')
	TransToGene('full.reads', 'full_gene.reads')