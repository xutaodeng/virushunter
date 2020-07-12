#!/usr/bin/env python
from collections import defaultdict
import sys
import string
import linecache

#consolidate two fq files
infile=sys.argv[1]
outfile1=sys.argv[2]
outfile2=sys.argv[3]
f1=open(infile, 'r')#'myseq2.fq', 'r')
#f2=open('f2.fq', 'r')
of1=open(outfile1, 'w')
of2=open(outfile2, 'w')


def Print(pair, lineno):
	global infile, of1, of2
	if pair=='1':
		print >> of1, linecache.getline(infile, lineno),
		print >> of1, linecache.getline(infile, lineno+1),
		print >> of1, linecache.getline(infile, lineno+2),
		print >> of1, linecache.getline(infile, lineno+3),
	else:
		print >> of2, linecache.getline(infile, lineno),
		print >> of2, linecache.getline(infile, lineno+1),
		print >> of2, linecache.getline(infile, lineno+2),
		print >> of2, linecache.getline(infile, lineno+3),

ids=defaultdict(list)
i=0
for line in f1:
	i+=1
	if i%4==1:
		id, pair=line.strip().split()
		pair=pair.split(':')[0]
		ids[id].append(pair+':'+str(i))
f1.close()

for id in ids.keys():
	if len(ids[id])==2:
		parts = ids[id]
		pair, lineno = parts[0].split(':')
		Print(pair, int(lineno))
		pair, lineno = parts[1].split(':')
		Print(pair, int(lineno))
of1.close()
of2.close()

# f1=open(infile, 'r')
# i=0
# for line in f1:
	# i+=1
	# if i%100000==0: print i
	# if i%4==1:
		# id, pair=line.strip().split()
		# pair=pair.split(':')[0]
		# if ids[id]==2:
			# doPrint=True
			# Print (line, pair)
		# else:
			# doPrint=False
	# elif doPrint:
		# Print (line, pair)
		
	
# f1.close()
		

