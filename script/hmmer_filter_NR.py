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
	f=open('/mnt/cluster/xdeng/blastdb/virus/virus.fa', 'r')
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
def readNRXML(fname, virusGIs):
	#print '---------------------------------------------------------------------'
	result_handle = open(fname, 'r')
	blast_records = NCBIXML.parse(result_handle)
	#print '---------------------------------------------------------------------'
	nrE={} # lowest E-values of non virus hit
	count=0 #count of nr in virusGIs
	for blast_record in blast_records:
		query = blast_record.query
		for alignment in blast_record.alignments:
			subject = alignment.title.upper()
			matches=re.findall(r'GI\|\d+\|', subject)
			v=False #virus is false
			for m in matches:
				if m in virusGIs:
					v=True
					break
			if v: count+=1; continue
			if ('VIRUS' in subject) or ('VIRAL' in subject): continue
			for hsp in alignment.hsps:
				if not nrE.has_key(query):
					nrE[query]=float(hsp.expect)
				elif float(hsp.expect) < float(nrE[query]):
					nrE[query]=float(hsp.expect)
	result_handle.close()
	#sys.stderr.write( 'count'+str( count))
	return nrE

# def readVirusXML1(fname):
	# result_handle = open(fname, 'r')
	# blast_records = NCBIXML.parse(result_handle)
	# virusE={}
	# virusID={}
	# for blast_record in blast_records:
		# query = blast_record.query
		# for alignment in blast_record.alignments:
			# subject = alignment.title.upper()
			# for hsp in alignment.hsps:
				# if not virusE.has_key(query):
					# virusE[query]=float(hsp.expect)
					# virusID[query]=subject # lowest subject
				# elif float(hsp.expect) < virusE[query]:
					# virusE[query]=float(hsp.expect)
					# virusID[query]=subject # lowest subject
	# result_handle.close()
	# return virusE, virusID

def hmmer_parser(fname, cachename):
	f=open(fname, 'r')
	start=False
	query=''
	subject=''
	virusE={}
	line=f.readline()
	while line:
		if line.strip().startswith ("Query:"):
			parts = line.strip().split()
			queryprot=parts[1]
			query = queryprot.rsplit('_', 1)[0][:-4]
			start=False
			end=False
		if not start and line.strip().startswith(">>"):
			start=True
			subject = line.strip()[3:].upper()
		elif line.strip().startswith("//") or line.strip().startswith(">>"):
			end=True
		if start and not end:
			if line.strip().startswith('== domain'):
				e_value = float(line.strip().split()[-1])
				score = int(float(line.strip().split()[-5]))
				a1=f.readline()
				a2=f.readline()
				a3=f.readline()
				if not virusE.has_key(query):
					query_nt = getSeq(cachename, cache, query)
					virusE[query]= [e_value, score, subject, a1,a2,a3,query_nt]
				elif e_value < virusE[query][0]:
					virusE[query][0:6]= e_value, score, subject,a1,a2,a3
		line=f.readline()
	f.close()
	return virusE

def hmmerOutputVirus(fname, filtertxt, hsp_only, E_VALUE_THRESH):
	global virusE, nrE, cache, cachename
	of = open(filtertxt, 'w')
	nalign, filter=0, 0
	for query in virusE.keys():
		if nrE.has_key(query) and virusE.has_key(query) and virusE[query][0] >= nrE[query]: 
			filter+=1; continue # filter out
		nalign+=1
		e_value, score, subject, a1,a2,a3, query_nt = virusE[query]
		of.write ('****Alignment****\n')
		of.write (query+'\n')
		query_nt = getSeq(cachename, cache, query)
		# if hsp_only=='YES':
			# ss,tt = int(hsp.query_start)-1, int(hsp.query_end)
			# query_nt = query_nt[ss:tt]
		of.write ( 'query_nt '+ query_nt+'\n')
		of.write ( 'subject: '+ subject+'\n')
		of.write ( 'length: '+ str(score)+'\n')
		of.write ( 'e value: '+ str(e_value)+'\n')
		try: nre= nrE[query]
		except: nre='no-hit'
		of.write ( 'lowest non-virus nr e value (LNVNRE) '+ str(nre)+'\n')
		of.write ( 'identities: '+ str(score)+'\n')
		of.write ( a1)
		of.write ( a2)
		of.write ( a3)
		# of.write ( str(hsp.query_start).ljust(11)+' '+ hsp.query+'\n')
		# of.write ( ' '.ljust(11)+' '+ hsp.match+'\n')
		# of.write ( str(hsp.sbjct_start).ljust(11)+' '+ hsp.sbjct+'\n')
	of.close()
	print fname, ' n_hits = ', str(nalign)
	print fname, ' nr_filter = ', str(filter) # filtered by nr

if __name__ == '__main__': 
	virusxml=sys.argv[1]
	nrxml=sys.argv[2]
	cachename = sys.argv[3]
	filtertxt = sys.argv[4]
	E_VALUE_THRESH = sys.argv[5]
	try: hsp_only = sys.argv[6]
	except: hsp_only = 'NO'
	#print 'hsp_only', hsp_only
	cache={}
	cache = CacheLines(cachename)
	#sys.stderr.write(cachename+'\n')
	virusGIs = readVirusGI()
	#sys.stderr.write('len(cache)'+str(len(cache)))
	virusE = hmmer_parser(virusxml, cachename)
	try: nrE = readNRXML(nrxml, virusGIs)
	except: nrE={}
	hmmerOutputVirus(virusxml, filtertxt, hsp_only, E_VALUE_THRESH)