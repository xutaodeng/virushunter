#!/usr/bin/env python

#this is the main 454 pipeline
from collections import defaultdict
import operator
import string
import sys
import os

#read both seqfile and qual file
def processfa(seqFile, qualFile, seqOut, qseqOut):
	print 'processing fa file', seqFile
	f = open(seqFile, 'r')
	qf = open(qualFile, 'r')
	of1=open(seqOut, 'w')
	of2=open(qseqOut, 'w')
	readID, seq, qseq=None, [], []
	i=0
	for line in f:
		qline=qf.readline()
		if line.strip().startswith('>'):
			i+=1
			if i%10000==0: print i
			if readID != None:
				of1.write('>'+readID+'\n'+''.join(seq)[26:]+'\n')
				qseqsub=qseq[26:]
				of2.write('>'+readID+'\n'+' '.join(qseqsub)+'\n')
			readID = line.strip().strip('>').split()[0]
			seq=[]
			qseq=[]
		else:
			seq.append(line.strip())
			qseq.extend(qline.strip().split())
	of1.write('>'+readID+'\n'+''.join(seq)[26:]+'\n')
	qseqsub=qseq[26:]
	of2.write('>'+readID+'\n'+' '.join(qseqsub)+'\n')
	f.close()
	qf.close()
	of1.close()
	of2.close()

if __name__ == "__main__":
	seqFile = 'nepal_Erie_in.sanger.fasta2' #'1.TCA.454Reads.fna'
	qualFile = 'nepal_Erie_in.sanger.fasta.qual2' #sys.argv[2] #'1.TCA.454Reads.qual'
	seqOut =  'nepal_Erie_in.sanger.fasta'  #sys.argv[3] #'1.TCA.454Reads.fna'
	qseqOut = 'nepal_Erie_in.sanger.fasta.qual' #sys.argv[4] #'1.TCA.454Reads.qual'
	processfa(seqFile, qualFile, seqOut, qseqOut)
