#!/usr/bin/env python

from collections import defaultdict
import operator
import sys
import os
import os.path

def PrintLines(fname, seqset, fsigname): 
	sys.stderr.write( 'Print...')
	f = open(fname, 'r')
	of = open(fsigname, 'w')
	i=0
	started=False
	for line in f:
		
		if i%10000000==0: sys.stderr.write(str(i)+'\n')
		if line.strip().startswith('>'):
			i+=1
			header=line.strip()[1:].strip()
			if started: started=False
			if header in seqset: 
				print >>of, line.strip() #header
				started=True
		elif started:
			print >>of, line.strip()
	of.close()
	f.close()

def screenBlast(fname): # extract subject or query titles 
	global LOW_E_VALUE_THRESH, HIGH_E_VALUE_THRESH
	sys.stderr.write('first screen blast\n')
	f = open(fname, 'r')
	seqset=set()
	for line in f:
		if line.startswith('#'): continue
		parts=line.strip().split()
		header=parts[1] #gnl|blast ID, index
		e = float(parts[-2])
		iden, alen = float(parts[2]), int(parts[3])
		if alen > 60 and HIGH_E_VALUE_THRESH>= e >= LOW_E_VALUE_THRESH:
			seqset.add(header)
	f.close()
	sys.stderr.write(str(len(seqset))+'\n')
	return seqset

if __name__ == '__main__': 
	fname=sys.argv[1] #xml file
	faname = sys.argv[2] # fa file
	fsigname=sys.argv[3]
	try: cachetype=sys.argv[4] #subject, default is query
	except: cachetype='query'
	HIGH_E_VALUE_THRESH=0.00001 #sys.argv[4]
	LOW_E_VALUE_THRESH=0 #sys.argv[4]
	seqset=screenBlast(fname)
	PrintLines(faname, seqset, fsigname)
