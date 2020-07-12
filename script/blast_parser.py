#!/usr/bin/env python

from Bio.Blast import NCBIXML
from collections import defaultdict
import operator
import sys
import os
import os.path
import linecache
import fcntl

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

if __name__ == '__main__': 
	fname=sys.argv[1]
	cachename = sys.argv[2]
	fsigname=sys.argv[3]
	mysfile=sys.argv[4]
	myslen=int(sys.argv[5])
	e_threshold=float(sys.argv[6])
	logfile=sys.argv[7]
	cache={}
	cache = CacheLines(cachename)
	result_handle = open(fname, 'r')
	
	queryset=set()
	sigReads=0
	f = open(fsigname, 'w')	
	
	try: 
		blast_records = NCBIXML.parse(result_handle)
		for blast_record in blast_records:
			for alignment in blast_record.alignments:
				for hsp in alignment.hsps:
					if hsp.expect < e_threshold:
						# print '****Alignment****'
						# print 'query', blast_record.query
						query_nt = getSeq(cachename, cache, blast_record.query)
						if blast_record.query not in queryset:
							if True:#(hsp.expect > 10E-15) :
								f.write('>'+blast_record.query+'\n')
								f.write(query_nt+'\n')
							sigReads+=1
						queryset.add(blast_record.query)
						# print 'query_nt', query_nt
						# print 'subject:', alignment.title
						# print 'length:', alignment.length
						# print 'e value:', hsp.expect
						# print 'identities:', hsp.identities
						# print str(hsp.query_start).ljust(11), hsp.query
						# print ' '.ljust(11), hsp.match
						# print str(hsp.sbjct_start).ljust(11), hsp.sbjct
	except:
		print 'bad xml', fsigname
		pass
	result_handle.close()
	f.close()
	nmys=print_mysterious(cachename, cache, queryset, mysfile, myslen)
	f10=open(logfile, 'a')
	fcntl.flock(f10, fcntl.LOCK_EX)
	f10.write(fname+' Sig_Vreads ='+' '+str(sigReads)+'\n')
	fcntl.flock(f10, fcntl.LOCK_UN)
	f10.close()