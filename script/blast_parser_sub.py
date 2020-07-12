#!/usr/bin/env python

from Bio.Blast import NCBIXML
from collections import defaultdict
import operator
import sys
import os
import os.path
import linecache

def CacheLines(fname, seqset): 
	sys.stderr.write( 'caching...')
	cache={}
	f = open(fname, 'r')
	i=0
	start, end = 0,0
	header= None
	for line in f:
		i+=1
		if i%10000000==0: sys.stderr.write(str(i))
		if line.strip().startswith('>'):
			end=i-1
			if header!=None and header in seqset: 
				cache[header] = (start, end)#;sys.stderr.write(header)
			header = line.strip()[1:].strip()
			start=i+1
	f.close()
	if header!=None and header in seqset: cache[header] = (start, i)
	sys.stderr.write('len(cache)'+str(len(cache))+'\n')
	return cache

def PrintLines(fname, seqset, fsigname): 
	sys.stderr.write( 'Print...')
	f = open(fname, 'r')
	of = open(fsigname, 'w')
	i=0
	started=False
	for line in f:
		i+=1
		if i%10000000==0: sys.stderr.write(str(i)+'\n')
		if line.strip().startswith('>'):
			if started: started=False
			header = line.strip()[1:].strip()
			if header in seqset: 
				print >>of, '>',header
				#sys.stderr.write(header+'\n')
				started=True
		elif started:
			print >>of, line.strip()
			#sys.stderr.write(line)
	of.close()
	f.close()

def getSeq(fsigname, cache, header):
	seq=[]
	start, end = cache[header]
	for i in xrange(start, end+1):
		seq.append(linecache.getline(fsigname, i).strip())
	return ''.join(seq)

def screenBlast(fname): # extract subject or query titles 
	global LOW_E_VALUE_THRESH, HIGH_E_VALUE_THRESH
	sys.stderr.write('first screen blast\n')
	result_handle = open(fname, 'r')
	blast_records = NCBIXML.parse(result_handle)
	seqset=set()
	f = open(fsigname, 'w')
	for blast_record in blast_records:
		for alignment in blast_record.alignments:
			subject=alignment.title.split(' ',1)[1].strip()
			for hsp in alignment.hsps:
				if HIGH_E_VALUE_THRESH>= hsp.expect >= LOW_E_VALUE_THRESH:
					seqset.add(subject)
					#print subject
					#sys.stderr.write( subject+'\n')
	result_handle.close()
	sys.stderr.write(str(len(seqset))+'\n')
	return seqset

if __name__ == '__main__': 
	fname=sys.argv[1] #xml file
	faname = sys.argv[2] # fa file
	fsigname=sys.argv[3]
	try: cachetype=sys.argv[4] #subject, default is query
	except: cachetype='query'
	HIGH_E_VALUE_THRESH=0.00001 #sys.argv[4]
	LOW_E_VALUE_THRESH=0 #sys.argv[4]
	seqset=screenBlast(fname)
	PrintLines(faname, seqset, fsigname)
	cache = CacheLines(fsigname, seqset)
	result_handle = open(fname, 'r')
	blast_records = NCBIXML.parse(result_handle)

	for blast_record in blast_records:
		for alignment in blast_record.alignments:
			subject=alignment.title.split(' ',1)[1].strip()
			for hsp in alignment.hsps:
				if  HIGH_E_VALUE_THRESH>= hsp.expect >= LOW_E_VALUE_THRESH:
					print '****Alignment****'
					if cachetype=='query':
						cache_seq = getSeq(fsigname, cache, blast_record.query)
					elif cachetype=='subject':
						cache_seq = getSeq(fsigname, cache, subject)
					print 'sequence', cache_seq
					print 'query', blast_record.query
					print 'subject:', alignment.title
					print 'length:', alignment.length
					print 'e value:', hsp.expect
					print 'identities:', hsp.identities
					print str(hsp.query_start).ljust(11), hsp.query
					print ' '.ljust(11), hsp.match
					print str(hsp.sbjct_start).ljust(11), hsp.sbjct
	result_handle.close()