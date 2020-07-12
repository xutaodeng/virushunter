#!/usr/bin/env python

from collections import defaultdict
import operator
import sys
import os
import os.path

f = open(sys.argv[1], 'r')
of = open(sys.argv[2], 'w')

line=f.readline().rstrip()
while line:
	if line.strip()=='<Iteration_hits>':
		of.write(line+'\n')
		line2=f.readline().strip()
		line3=f.readline().strip()
		if not (line2=='</Hit_hsps>' and line3=='</Hit>'): 
			of.write(line2+'\n')
			of.write(line3+'\n')
	else:
		of.write(line+'\n')
	line=f.readline().rstrip()
print 'done'
f.close()
of.close()
