#!/usr/bin/env python

from collections import defaultdict
import operator
import sys
import os
import re

def trimf2ends(infile, outfile): #illumina 33 or 64
	f = open(infile, 'r')
	of = open(outfile, 'w')
	i = 0
	n5,n3, goodread=0,0,0
	totalquality, trimquality=0,0
	totallength, trimlength=0,0
	for line in f:
		i+=1
		#if i==100000: break
		if i%4==1:
			id=line.strip()
		elif i%4==2: #sequence string
			read=line.strip()
		elif i%4==3:
			qid=line.strip()
		elif i%4==0: #quality string
			qseq=line.strip()
			newread=[]
			phreds=[ord(q)-33 for q in qseq]
			totalquality+=sum(phreds)
			totallength+=len(read)
			newread=''.join([r if p >= 30 else 'N' for p, r in zip(phreds, read)])
			s,t=0,0
			try: s = len(re.findall('\AN+', newread)[0])
			except: pass
			try: t = len(re.findall('N+$', newread)[0])
			except: pass 
			if s>0: n5+=1
			if t>0: n3+=1
			t=len(newread)-t
			nr=newread[s:t]
			nq=qseq[s:t]
			trimquality+=sum([ord(q)-33 for q in nq])
			trimlength+=len(nq)
			if len(nr)>30 and nr.count('N')/float(len(nr))<0.1:
				goodread+=1
				print >>of, '\n'.join([id, nr, qid, nq])
	f.close()
	of.close()
	print infile
	print 'number of 5 primers', n5
	print 'number of 3 primers', n3
	print 'number of reads', i/4, 'good reads', goodread, 'good rate', goodread*4/float(i)*100 
	print 'average quality before trim', totalquality/totallength, 'after trim', trimquality/trimlength

if __name__ == '__main__':
	infile = sys.argv[1]
	outfile =sys.argv[2] 
	trimf2ends(infile, outfile)

