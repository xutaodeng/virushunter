#!/usr/bin/env python

#this is the main 454 pipeline
from collections import defaultdict
import operator
import string
import sys
from Bio import pairwise2
from multiprocessing import Pool
import os

def revcomp(tag):
	rval = []
	for b in tag:
		if b == 'A': a='T'
		elif b == 'C': a='G'
		elif b == 'G': a='C'
		elif b == 'T': a='A'
		try: rval.append(a)
		except: print 'error', tag
	rval.reverse()
	return ''.join(rval)
	
def openBarcodes(tagFile):
	global wd, base
	ff=open(tagFile, 'r')
	tags={}
	fhandles={}
	header = ff.readline()
	for line in ff:
		parts = line.strip().split()
		id, tag = parts[0], parts[1]
		print 'tagid', id
		tags[id]= (tag, revcomp(tag))
		fhandles[id]=open(id+'_.fastq', 'w')
	fhandles['undetermined']=open('undetermined_.fastq', 'w')
	ff.close()
	return tags, fhandles

def closeBarcodes(fhandles):
	for (id, handle) in fhandles.items():
		handle.close()

#read both seqfile and qual file
def processfq(seqFile):
	print 'processing fq file', seqFile
	f = open(seqFile, 'r')
	i=0
	while 1:
		i+=1
		id=f.readline()
		if i%1000==1: print id
		if not id: break
		seq=f.readline()
		qid=f.readline()
		qseq=f.readline()
		readID= processRead(seq)
		fhandles[readID].write(id+seq[20:]+'\n'+qid+qseq[20:]+'\n')
	f.close()

def processRead(seq):
	global tags
	rval='undetermined'
	for (id, tag) in tags.items():
		ftag, rtag = tag[0], tag[1]
		fseq=seq[0:20]
		rseq=seq[-20:]
		score1 = pairwise2.align.localms(fseq, ftag, 1, 0, -5, -1, score_only=True, one_alignment_only=True)
		score2 = pairwise2.align.localms(rseq, rtag, 1, 0, -5, -1, score_only=True, one_alignment_only=True)
		if score1>=18:
			rval =id
		elif score2>=18:
			rval =id
	return rval

# def generateMira(seqFile):
	# global tags, wd, base
	# of =open(os.path.abspath(wd)+'/mira.sh', 'w')
	# for id in tags.keys():
		# of.write('mira --project='+base+'_'+id+' --job=denovo,genome,accurate,sanger -MI:sonfs=no --fasta\n')
	# of.close()

if __name__ == "__main__":

	seqFile = sys.argv[1] 
	tagFile= sys.argv[2] #sys.argv[3] #'Colorado.Res.Infect..txt'

	tags, fhandles = openBarcodes(tagFile)
	processfq(seqFile)
	closeBarcodes(fhandles)


