#!/usr/bin/env python
from collections import defaultdict
import os.path
import os
import string
import sys
from os import walk
from os import listdir

def combine_result(labels):
	fs=[]
	of=open('gene_ex_combine_diff.txt', 'w')
	for label in labels:
		f=open(label+'/diff_d_out/gene_exp.diff', 'r')
		fs.append(f)
	j=0
	for line in fs[0]:
		if j==0: #header
			line=label+'_'+line
		lines=[line.strip()]
		for i in range(1, len(fs)):
			li=fs[i].readline().strip()
			if j==0: #header
				li=label+'_'+li
			lines.append(li)
		out='\t'.join(lines)
		of.write(out+'\n')
		j+=1

	for f in fs:
		f.close()
	of.close()

if __name__ == "__main__":
	labels = sys.argv[1:]
	combine_result(labels)