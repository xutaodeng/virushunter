#!/usr/bin/env python
from Bio.Blast import NCBIXML
from collections import defaultdict
import operator
import sys
import os
import os.path
import linecache
import re

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

def readVirusGI():
	f=open('/mnt/cluster/xdeng/blastdb/virus.fa', 'r')
	virusGIs=[]
	for line in f:
		if line.strip().startswith('>'):
			try: 
				gi = line.strip().split('|')[1]
				virusGIs.append('GI|'+gi+'|')
			except:
				pass
	return set(virusGIs)

def getSeq(cachename, cache, header):
	seq=[]
	start, end = cache[header]
	for i in xrange(start, end+1):
		seq.append(linecache.getline(cachename, i).strip())
	return ''.join(seq)

def readDiamondNR(fname):
	f = open(fname, 'r')
	nrE=set([]) # nr blacklist
	for line in f:
		parts = line.strip().split()
		query, subject = parts[0], parts[1]
		subject=subject.strip().split('_', 1)[0]
		if subject != "VIRUS":
			nrE.add(query)
	f.close()

	return nrE

def readVirusXML1(fname):
	result_handle = open(fname, 'r')
	virusE={}
	virusID={}
	try:
		blast_records = NCBIXML.parse(result_handle)
		for blast_record in blast_records:
			query = blast_record.query
			for alignment in blast_record.alignments:
				subject = alignment.title.upper()
				for hsp in alignment.hsps:
					if not virusE.has_key(query):
						virusE[query]=float(hsp.expect)
						virusID[query]=subject # lowest subject
					elif float(hsp.expect) < virusE[query]:
						virusE[query]=float(hsp.expect)
						virusID[query]=subject # lowest subject
	except:
		pass
	result_handle.close()
	return virusE, virusID

def OutputVirus(fname, filtertxt, hsp_only, E_VALUE_THRESH):
	global virusE, nrE, cache, cachename
	of = open(filtertxt, 'w')
	#library=os.path.basename(filtertxt).rsplit('_',3)[0] #'Ita-Res-C11-15_blast_filter.txt_34'
	result_handle = open(fname, 'r')
	blast_records = NCBIXML.parse(result_handle)
	nalign, filter=0, 0
	try: 
		for blast_record in blast_records:
			#E_VALUE_THRESH = 0.01
			query = blast_record.query
			for alignment in blast_record.alignments:
				for hsp in alignment.hsps:
					if float(hsp.expect) < E_VALUE_THRESH:
						#try: 
							#sys.stderr.write('nrE: '+str(nrE[query])+' vrE: '+str(virusE[query])+'\n')
							#sys.stderr.write('nrID: '+str(nrID[query])+'\n')
							#sys.stderr.write('virus: '+alignment.title+'\n')
						#except: pass
						if query in nrE: 
							filter+=1; continue # filter out
						nalign+=1
						of.write ('****Alignment****\n')
						of.write ( blast_record.query+'\n')
						query_nt = getSeq(cachename, cache, blast_record.query)
						if hsp_only=='YES':
							ss,tt = int(hsp.query_start)-1, int(hsp.query_end)
							query_nt = query_nt[ss:tt]
							#print ss, tt, query_nt
						of.write ( 'query_nt '+ query_nt+'\n')
						of.write ( 'subject: '+ alignment.title+'\n')
						of.write ( 'length: '+ str(alignment.length)+'\n')
						of.write ( 'e value: '+ str(hsp.expect)+'\n')
						# try: nre= nrE[query]
						# except: nre='no-hit'
						of.write ( 'lowest non-virus nr e value (LNVNRE) '+ str('-')+'\n')
						of.write ( 'identities: '+ str(hsp.identities)+'\n')
						of.write ( str(hsp.query_start).ljust(11)+' '+ hsp.query+'\n')
						of.write ( ' '.ljust(11)+' '+ hsp.match+'\n')
						of.write ( str(hsp.sbjct_start).ljust(11)+' '+ hsp.sbjct+'\n')
	except: print 'XML format bad' 
	result_handle.close()
	of.close()
	print fname, ' n_hits = ', str(nalign)
	#print fname, ' n_filter = ', str(filter)

if __name__ == '__main__': 
	virusxml=sys.argv[1]
	nrDiamond=sys.argv[2]
	cachename = sys.argv[3]
	filtertxt = sys.argv[4]
	E_VALUE_THRESH = sys.argv[5]
	try: hsp_only = sys.argv[6]
	except: hsp_only = 'NO'
	print 'hsp_only', hsp_only
	cache={}
	cache = CacheLines(cachename)
	#sys.stderr.write(cachename+'\n')
	virusGIs = readVirusGI()
	#sys.stderr.write('len(cache)'+str(len(cache)))
	virusE, virusID = readVirusXML1(virusxml)
	try: nrE = readDiamondNR(nrDiamond)
	except: nrE=set([]); print 'no diamond'
	OutputVirus(virusxml, filtertxt, hsp_only, E_VALUE_THRESH)
