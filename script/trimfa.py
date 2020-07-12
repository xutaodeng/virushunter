#!/usr/bin/env python
from collections import defaultdict
import operator
import sys
import os
import os.path

def filterLength(infile, outfile, length):
	print 'converting', infile, 'to', outfile
	f = open(infile, 'r')
	of = open(outfile, 'w')
	i = 0
	seq=[]
	id=''
	a, b = 0, 0
	for line in f:
		if line.strip().startswith('>'):
			a+=1
			if id!='' and len (''.join(seq)) >= length:
				b+=1
				of.write('>'+id+'\n'+''.join(seq)+'\n')
			id = line.strip()[1:]
			seq=[]
		else:
			seq.append(line.strip())
	of.write('>'+id+'\n'+''.join(seq)+'\n')
	f.close()
	of.close()

if __name__ == '__main__':
	infile, outfile, length =sys.argv[1], sys.argv[2], sys.argv[3]
	#outfile = infile.rsplit('.', 1)[0]+str(l)+'_filter.fa'
	filterLength(infile, outfile, int(length))