#!/usr/bin/env python
from collections import defaultdict
import operator
import sys
import os
import os.path
import itertools

def removedup(inputfq, outputfq):
	#print 'processing..', inputfq
	f=open(inputfq, 'r')
	of=open(outputfq, 'a')
	i=0
	keyset=set()
	for line in f:
		i+=1
		if i%4==1:
			id = line.strip()
		elif i%4==2:
			read = line.strip()
		elif i%4==3:
			id2 = line.strip()
		elif i%4==0:
			qual = line.strip()
			if read[5:45] in keyset: pass
			else:
				print >>of, '\n'.join([id,read,id2,qual])
				keyset.add(read[5:45])
	f.close()
	of.close()

def split(filename):
	print 'spliting..'
	comb=itertools.product(['A','C','G','T'], repeat=4)
	keydict={}
	for it in comb: 
		key=''.join(it)
		keydict[key]=open(filename+'_'+key,'w')

	i=0
	for line in f:
		i+=1
		if i%4==1:
			id = line.strip()
		elif i%4==2:
			read = line.strip()
			readstart=read[0:4]
		elif i%4==3:
			id2 = line.strip()
		elif i%4==0:
			qual = line.strip()
			try: print >>keydict[readstart], '\n'.join([id,read,id2,qual])
			except: pass
	f.close()
	for it in comb:
		key=''.join(it)
		keydict[key].close()
	of.close()
def cleanTempFiles():
	comb=itertools.product(['A','C','G','T'], repeat=4)
	keydict={}
	for it in comb: 
		key=''.join(it)
		os.remove(filename+'_'+key)
if __name__ == '__main__':
	try:
		filename=sys.argv[1]
		f=open(filename, 'r')
		outfile=sys.argv[2]
		of=open(outfile, 'w') #clear this output file
		of.close()
		split(filename)

		print 'clonal reads removing...'
		comb=itertools.product(['A','C','G','T'], repeat=4)
		for it in comb: 
			key=''.join(it)
			removedup(filename+'_'+key, outfile)
		cleanTempFiles()
	except (KeyboardInterrupt, SystemExit):
		cleanTempFiles()

	