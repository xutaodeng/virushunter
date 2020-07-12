#!/usr/bin/env python
from collections import defaultdict
import operator
import sys
import os
import os.path
import itertools
import gzip
import random

def polyA(filename, keypath, raw):
	readlen=0
	if filename.endswith('.gz'):
		f=gzip.open(filename, 'r')
	else:
		f=open(filename, 'r')
	hist = keypath+'_'+raw+'_hist.txt'
	histpng=keypath+'_'+raw+'_hist.png'
	histR=keypath+'_'+raw+'_hist.R'
	outfile=keypath+'_'+raw+'_hist.Rout'
	key = os.path.basename(keypath)
	title=key+'_'+raw
	f2=open(hist, 'w')
	polya, polyt = 0, 0
	i=0
	totalRead, uniRead=0, 0
	maxRead=0
	length=[]
	for line in f:
		i+=1
		if i%4==2:
			totalRead+=1
			read = line.strip()
			if len(read)>=20:
				if len(read)> maxRead: maxRead=len(read)
				uniRead+=1
				readlen+=len(read)
				length.append(len(read))
				if 'AAAAA' in read: polya+=1
	f.close()
	try: sel = random.sample(length, 1000)
	except: sel = random.sample(length, len(length))
	for ss in sel:
		f2.write(str(ss)+'\n')
	f2.close()
	if raw!='raw':
		print filename, 'uni_reads =', uniRead
		print filename, 'clean_read_length =', float(readlen)/uniRead
		print filename, 'num_polyA_reads =', polya
	else:
		print filename,  'total_reads = ', totalRead
		print filename, 'raw_read_length =', maxRead
	of3=open(histR, 'w')
	Rscript =' a<-read.table("'+hist+'", header=F,  stringsAsFactors=F)\n'+ \
			'family <- as.numeric(a[[1]])\n' +\
			'png("'+histpng+'", res=90, width=600, height=600)\n' +\
			'maxlen<-max(family)+20\n' +\
			'hist(family, xlab="read length", main="'+title+'", col="red",  breaks=seq(20, maxlen, 5))\n'+\
			'dev.off()\n'
	of3.write(Rscript)
	of3.close()
	cmd='R CMD BATCH --quiet --vanilla '+histR+' '+outfile
	print cmd
	os.system(cmd)

if __name__ == '__main__':
	filename = sys.argv[1]
	keypath = sys.argv[2]
	try: raw=sys.argv[3]
	except: raw='clean'
	polyA(filename, keypath, raw)