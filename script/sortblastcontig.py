#!/usr/bin/env python
import sys
from collections import defaultdict
from operator import itemgetter



infile =sys.argv[1]
f =open(infile, 'r')
dict=defaultdict(list)
for line in f:
	if line.strip().startswith('blastn'): continue
	else: 
		program, key, vid, cov = line.strip().split()
		gi=vid.split('|')[1]
		dict[key+':'+gi].append((program, int(cov)))
f.close()
for (kk, value) in dict.items():
	v =sorted(value, key=itemgetter(1), reverse=True)
	if v[0][1] > 500: 
		for dd in v: print kk, dd[0], dd[1]