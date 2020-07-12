#!/usr/bin/env python

from collections import defaultdict
import operator
import sys
import os
import os.path
import linecache

def CacheLines(fname): 
	cache={}
	f = open(fname, 'r')
	i=0
	start, end = 0,0
	header= None
	for line in f:
		i+=1
		if line.strip().startswith('>'):
			end=i-1
			if header!=None: cache[header] = (start, end)
			header = line.strip()[1:]
			start=i+1
	if header!=None: cache[header] = (start, i)
	return cache

def getSeq(cachename, cache, header):
	seq=[]
	start, end = cache[header]
	for i in xrange(start, end+1):
		seq.append(linecache.getline(cachename, i).strip())
	return ''.join(seq)

def print_mysterious(cachename, cache, queryset, outfile, length):
	nmys=0
	of=open(outfile, 'w')
	for header in cache.keys():
		if header in queryset: #hits
			continue
		seq=getSeq(cachename, cache, header)
		if len(seq)> length:
			nmys+=1
			of.write('>'+header+'\n')
			of.write(seq+'\n')
	of.close()
	return nmys #number of mysterious contigs

def hmmer_parser(fname, cachename, fsigname):
	f=open(fname, 'r')
	of = open(fsigname, 'w')
	start=False
	queryset=set()
	querydna=''
	sigReads=0
	for line in f:
		if line.strip().startswith ("Query:"):
			if start and querydna not in queryset: 
				query_nt = getSeq(cachename, cache, querydna)
				of.write('>'+querydna+'\n')
				of.write(query_nt+'\n')
				sigReads+=1
				queryset.add(querydna)
			parts = line.strip().split()
			queryprot=parts[1]
			querydna = queryprot.rsplit('_', 1)[0][:-4]
			start=False
			# end=False
			# record=[]
		if not start and line.strip().startswith(">>"):
			start=True
		# elif line.strip().startswith("//") or line.strip().startswith(">>"):
			# end=True
		# if start and not end:
			# record.append(line.strip())
	if start and querydna not in queryset: 
		query_nt = getSeq(cachename, cache, querydna)
		of.write('>'+querydna+'\n')
		of.write(query_nt+'\n')
		sigReads+=1
	f.close()
	of.close()
	print fname, 'Sig_Vreads =', sigReads
	return queryset

if __name__ == '__main__': 
	fname=sys.argv[1]
	cachename = sys.argv[2]
	fsigname=sys.argv[3]
	mysfile=sys.argv[4]
	myslen=int(sys.argv[5])
	e_threshold=float(sys.argv[6])
	cache = CacheLines(cachename)
	queryset = hmmer_parser(fname, cachename,fsigname)
	nmys=print_mysterious(cachename, cache, queryset, mysfile, myslen)