#!/usr/bin/env python
import sys
import os
from Bio import pairwise2
def GetN50(lenc):
	lenc=sorted(lenc, reverse=True)
	total=sum(lenc)
	cum=0
	for le in lenc:
		cum+=le
		if cum*2 > total:
			return le
def Similarity(s1, s2):
	score = pairwise2.align.localms(s1, s2, 2, -1, -2, -1, score_only=True)
	return score, len(s1), len(s2)

def readContig(filename, len_thres, ofile):
	diclen={}
	f=open(filename, 'r')
	of=open(ofile, 'w')
	id=''
	contig=[]
	seqs={}
	s=0
	m=0 # max
	lenc=[]
	i=0
	for line in f:
		if line.strip().startswith('>'):
			contig=''.join(contig)
			if len(contig) >len_thres:
				i+=1
				seqs[str(i)]=''.join(contig)
				lenc.append(len(contig))
				of.write('>'+str(i)+'\n'+contig+'\n')
			id=line.strip()
			contig=[]
		else:
			contig.append(line.strip())
	if len(contig) >len_thres:
		seqs[str(i)]=''.join(contig)
		lenc.append(len(contig))
		of.write(line.strip()+'\n'+contig+'\n')
	f.close()
	of.close()
	ncontigs=len(seqs)
	total=sum(lenc)
	ma=max(lenc)
	ave = total/ncontigs
	N50 = GetN50(lenc)
	return ncontigs, ma, ave, N50, seqs

file1=sys.argv[1]
file2=sys.argv[2]
ncontigs1, ma1, ave1, N501, seqs1 = readContig(file1, 1000, 'test1.fa')
ncontigs2, ma2, ave2, N502, seqs2 = readContig(file2, 1000, 'test2.fa')

print 'file1', ncontigs1, ma1, ave1, N501
print 'file2', ncontigs2, ma2, ave2, N502

cmd='blastn -task megablast -max_target_seqs 1 -outfmt 6 -out test.tab -query test1.fa -subject test2.fa'
print cmd
os.system(cmd)
f=open('test.tab', 'r')
match1=0
diclen1, diclen2=set(), set()
for line in f:
	parts = line.strip().split()
	id1, id2 = parts[0], parts[1]
	diclen1.add(id1)
	diclen2.add(id2)
	le = int(parts[3])
	match1+=le
f.close()
total1, total2= 0, 0
for key in seqs1.keys():
	if key in diclen1: pass
	total1 += len(seqs1[key])
for key in seqs2.keys():
	if key in diclen2: pass
	total2 += len(seqs2[key])
print match1, total1, total2