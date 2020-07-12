#!/usr/bin/env python
from Bio.Blast import NCBIXML
from collections import defaultdict
import operator
import sys
import os
import os.path
import linecache
import re

def FilterLines(submyspre, submys, hits): 
	f = open(submyspre, 'r')
	of = open(submys, 'w')
	nmys=0
	out=False
	for line in f:
		if line.strip().startswith('>'):
			header = line.strip()[1:].strip()
			if header not in hits:
				of.write(line.strip()+'\n')
				nmys+=1
				out=True
		elif out:
			of.write(line.strip()+'\n')
			out=False
	#print submys, 'n_mys =', nmys

def ReadXML(submysxmlpre):
	hits=set([])
	try: 
		result_handle = open(submysxmlpre, 'r')
		blast_records = NCBIXML.parse(result_handle)
		for blast_record in blast_records:
			for alignment in blast_record.alignments:
				for hsp in alignment.hsps:
					hits.add(blast_record.query)
		result_handle.close()
	except: return hits
	return hits

if __name__ == '__main__': 
	submysxmlpre=sys.argv[1] #xml pre mys
	submyspre = sys.argv[2] #pre mys fasta
	submys = sys.argv[3] #mys fasta filtered
	hits = ReadXML(submysxmlpre)
	FilterLines(submyspre, submys, hits)