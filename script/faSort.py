#!/usr/bin/env python
import sys
from operator import itemgetter, attrgetter
import linecache

def filterSort(infile):
	length=[]
	f = open(infile, 'r')
	seq=[]
	id=''
	i=0
	start, end =0,0
	nlower=0
	for line in f:
		i+=1
		if line.strip().startswith('>'):
			end = i-1
			if id!='' : 
				seqlen = len (''.join(seq))
				length.append((start, end, seqlen, nlower))
			id = line.strip()
			seq=[]
			start=i
			end=0
			nlower=0
		else:
			ss=line.strip()
			nlower += len([c for c in ss if c.islower()])
			seq.append(ss)
	end = i
	seqlen = len (''.join(seq))
	if id!='': length.append((start,end, seqlen, nlower))
	f.close()
	
	length.sort(key=itemgetter(2), reverse=True)
	return length

def printFa(infile, outfile):
	of = open(outfile, 'w')
	nmys=0
	for (start, end, seqlen, nlower) in length:
		#if nlower/float(seqlen) > 0.2: continue
		nmys+=1
		for i in xrange(start, end+1):
			line = linecache.getline(infile, i)
			if i==start: of.write(line.strip()+' length='+str(seqlen)+'\n')
			else: of.write(line)
	of.close()
	if log=='True':
		print outfile, 'n_mys =', nmys

if __name__ == '__main__':
	infile, outfile= sys.argv[1], sys.argv[2]
	try: log=sys.argv[3]
	except: log='False'
	length = filterSort(infile)
	printFa(infile, outfile)
