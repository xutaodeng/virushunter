#!/usr/bin/env python

from collections import defaultdict
import operator
import sys
import os
import re

def readFa(input, output, mode):
	f=open(input, 'r')
	of = open(output, mode)
	rval={}
	dup=False
	for line in f:
		if line.strip().startswith('>'):
			id, annot=line.strip()[1:].split('dummy',1)
			if not rval.has_key(id):
				dup=False
				of.write('>'+annot+':'+id+'\n') #id
			else: dup =True
			rval[id]=annot
		elif not dup:
			of.write(line)
	return rval
	
def annotate(tally, hits, mystery, out, hitfa, mysteryfa): 
	f1 =open(tally, 'r')
	hitid=readFa(hits, hitfa, 'w')
	print 'hitid', len(hitid)
	mysid=readFa(mystery, mysteryfa, 'w')
	print 'mysid', len(mysid)
	of =open(out, 'w')
	nhit,nmys=0,0
	header=f1.readline()
	of.write('type\t'+header)
	for line in f1:
		id=line.strip().split()[0].split('dummy',1)[0]
		print id
		if hitid.has_key(id):
			of.write(hitid[id]+'\t'+line.strip().replace('_dummy','')+'\n')
			nhit+=1
		elif mysid.has_key(id):
			of.write('mysterious\t'+line.strip().replace('_dummy','')+'\n')
			nmys+=1
	f1.close()
	of.close()
	print 'nhit', nhit, 'nmys', nmys

if __name__ == '__main__':
	tally, hits, mystery, out, hitfa, mysteryfa = sys.argv[1:]
	annotate(tally, hits, mystery, out, hitfa, mysteryfa )
