#!/usr/bin/env python

from Bio.Blast import NCBIXML
from collections import defaultdict
import operator
import sys
import os
import os.path

if __name__ == '__main__': 
	fname=sys.argv[1] #'blast_out.xml'
	fsigname=sys.argv[2] #'blast_extract.fasta'
	batchfile=sys.argv[3] # 'batch.txt'
	
	result_handle = open(fname, 'r')
	blast_records = NCBIXML.parse(result_handle)

	queryset=set()
	sigReads=0
	f = open(fsigname, 'w')
	f2=open(batchfile, 'w')
	for blast_record in blast_records:
		#E_VALUE_THRESH = e_threshold
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				if True: #hsp.expect < E_VALUE_THRESH:
					# print '****Alignment****'
					# print 'query', blast_record.query
					# query_nt = getSeq(cachename, cache, blast_record.query)
					# if blast_record.query not in queryset:
						# f.write('>'+blast_record.query+'\n')
						# f.write(query_nt+'\n')
					sigReads+=1
					# queryset.add(blast_record.query)
					# print 'query_nt', query_nt
					# print 'subject:', alignment.title
					extract=alignment.title.strip().split()[0]
					if extract in queryset: continue
					queryset.add(extract)
					f.write('>'+alignment.title+'\n')
					f2.write(extract+'\n')
					f.write(hsp.sbjct.replace('-','')+'\n')
					# print 'length:', alignment.length
					# print 'e value:', hsp.expect
					# print 'identities:', hsp.identities
					# print str(hsp.query_start).ljust(11), hsp.query
					# print ' '.ljust(11), hsp.match
					# print str(hsp.sbjct_start).ljust(11), hsp.sbjct
	result_handle.close()
	f.close()
	f2.close()
	print 'num_alignment', sigReads
