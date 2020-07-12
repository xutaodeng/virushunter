#!/usr/bin/env python

from collections import defaultdict, deque
import operator
import sys
import os
import re

def complex(seq):
	if len(seq)<30: return True #too short
	words=deque()
	wordfreq=defaultdict(int)
	maxw=0
	maxword=set([])
	for i in xrange(len(seq)-2):
		word=seq[i:i+3]
		if word=='NNN': return True
		words.append(word)
		worfreq[word]+=1
		if wordfreq[word]>maxw: maxw+=1; maxword.add(word)
		if len(words)==20 and len(set(words)) <= 3:
			#print seq[(i-19):(i+3)]#, len(set(words)), words
			return True
		if len(words)>20:
			w=words.popleft()
			wordfreq[w]-=1
			if maxword
			
	return False

def process(infile, outfile): #illumina 33 or 64
	f = open(infile, 'r')
	of = open(outfile, 'w')
	i = 0
	good, bad=0,0
	for line in f:
		i+=1
		if i==10000: print i
		if i%4==1:
			id=line.strip()
		elif i%4==2: #sequence string
			read=line.strip()
		elif i%4==3:
			qid=line.strip()
		elif i%4==0: #quality string
			qseq=line.strip()
			if len(read)<30 : continue
			if complex(read): bad+=1; continue
			good+=1
			print >>of, '\n'.join([id, read, qid, qseq])
	f.close()
	of.close()
	print infile, 'good', good, 'bad', bad


if __name__ == '__main__':
	infile = sys.argv[1]
	outfile =sys.argv[2]
	process(infile, outfile)