#!/usr/bin/env python

from Bio.Blast import NCBIXML
from collections import defaultdict
import operator
import sys
import os
import os.path

if __name__ == '__main__': 
	fname=sys.argv[1]
	fsigname=sys.argv[2]
	f=open(fsigname, 'w')
	result_handle = open(fname, 'r')
	blast_records = NCBIXML.parse(result_handle)
	f = open(fsigname, 'w')
	for blast_record in blast_records:
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				f.write('>'+alignment.title+'\n')
				f.write(hsp.sbjct.replace('-', '')+'\n')
	result_handle.close()
	f.close()
