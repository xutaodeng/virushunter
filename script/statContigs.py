#!/usr/bin/env python
import sys
import os

def statContigs(lenc):
	lenc=sorted(lenc, reverse=True)
	total=sum(lenc)
	if len(lenc)>=3: top3= [lenc[0], lenc[1], lenc[2]]
	elif len(lenc)>=2: top3= [lenc[0], lenc[1], 0]
	elif len(lenc)>=1: top3= [lenc[0], 0, 0]
	else: top3= [0,0,0]
	cum=0
	N50=0
	n500, n1000, n2000,n3000, n4000,n5000,n10000=0, 0, 0, 0,0, 0,0
	for le in lenc:
		if le>500: n500+=le
		if le>1000: n1000+=le
		if le>2000: n2000+=le
		if le>3000: n3000+=le
		if le>4000: n4000+=le
		if le>5000: n5000+=le
		if le>10000: n10000+=le

	for le in lenc:
		cum+=le
		if cum*2 > total:
			N50=le
			break
	
	ncontigs=len(lenc)
	total=sum(lenc)
	# ma=max(lenc)
	# ave = total/ncontigs
	#N50 = GetN50(lenc)
	return top3, N50, ncontigs, n500/1000.0, n1000/1000.0, n2000/1000.0, n5000/1000.0

def GetN50(lenc):
	lenc=sorted(lenc, reverse=True)
	total=sum(lenc)
	cum=0
	for le in lenc:
		cum+=le
		if cum*2 > total:
			return le
	return 0

def readContig(filename, thres=200):
	diclen={}
	f=open(filename, 'r')
	contig=[]
	seqs={}
	s=0
	m=0 # max
	lenc=[]
	for line in f:
		if line.strip().startswith('>'):
			contig_length = len(''.join(contig))
			if contig_length >= thres: lenc.append(contig_length)
			contig=[]
		else:
			contig.append(line.strip())
	if contig!=[]:
		contig_length = len(''.join(contig))
		if contig_length >= thres: lenc.append(contig_length)
	f.close()
	return lenc


file1=sys.argv[1]
program = sys.argv[2]
key = sys.argv[3]
try: label =sys.argv[4]
except: label = 'none'
lenc=readContig(file1)
top3, N50, ncontigs, n500, n1000, n2000,n5000 = statContigs(lenc)

print program, label, key, top3[0], top3[1], top3[2], N50, ncontigs, n500, n1000, n2000, n5000
#print program key top1 top2 top3 N50 ncontigs n500 n1000 n2000 n5000 n10000