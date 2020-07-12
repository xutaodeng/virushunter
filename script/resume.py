#!/usr/bin/env python
from collections import defaultdict
import operator
import sys
import os
import os.path
import itertools
import gzip
import random

def resume(shfile):
	f=open(shfile, 'r')
	of =open('newscript.sh', 'w')
	i=0
	for line in f:
		if line.strip()=='wait':
			of.write(line.strip()+'\n')
			continue
		parts = line.strip().split()
		filename=''
		for part in parts:
			if 'xml' in part:
				filename= part
				break
		if os.path.isfile(filename) and os.stat(filename).st_size > 0:
			print filename
			continue
		of.write(line.strip()+'\n')
		i+=1
	f.close()
	of.close()
	print i, 'files to run'

if __name__ == '__main__':
	shfile = sys.argv[1]
	resume(shfile)