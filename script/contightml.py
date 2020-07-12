#!/usr/bin/env python

from collections import defaultdict
import operator
import sys
import os

if __name__ == '__main__':
	infile = sys.argv[1] #'s_combine_2_sequence.txt'
	outhtml =sys.argv[2] #'s_combine_trim_2_sequence.txt'
	f = open(infile, 'r')
	of = open(outhtml, 'w')
	for line in f:
		parts = line.strip().split(':', 1)
		if 'max_contig' in line:
			#contig_AaC/contig.blast:cap3_SOcaseDRNoAmp23 fecalphage8L3 max_contig:13809<br>
			out1 = parts[0].split('/')[0].split('_')[1]
			out2 = parts[1].split('_', 1)[1]
			out=out1+' '+out2
			if 'individual' in line:
				out = parts[1]
			of.write(out+'\n')
		else:
			of.write(parts[1]+'\n')
	f.close()
	of.close()
