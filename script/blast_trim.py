#!/usr/bin/env python
import sys
from collections import defaultdict
from operator import itemgetter, attrgetter

f=open(sys.argv[1], 'r')#fq file
f2=open(sys.argv[2],'r') #blast table
of=open(sys.argv[3]+'.tmp','w') #output fq file
label =sys.argv[4] #library key
try: pair_end =sys.argv[5] #pair 1 or 2
except: pair_end ='1'

def revcomp(seq):
	rval=[]
	for s in seq:
		if s=='A':
			rval.append('T')
		elif s=='T':
			rval.append('A')
		elif s=='G':
			rval.append('C')
		elif s=='C':
			rval.append('G')
		else:
			rval.append(s)
	rval.reverse()
	return ''.join(rval)

def trimAdaptorEnd(infile, outfile, ada): #input fq file and output fq file
	f=open(infile, 'r')
	of=open(outfile, 'w')
	rada=revcomp(ada)
	partial=0
	i=0
	for line in f:
		i+=1
		if i%4==2:
			seq=line.strip()
			rval1=checkEnd(seq, ada)
			rval2=checkEnd(seq, rada)
			rval=max(rval1, rval2) #rval -1 and above
			if rval>=0: partial+=1
			if len(seq[(rval+1):])==0: of.write('A\n')
			else: of.write(seq[(rval+1):]+'\n')
		elif i%4==0:
			if len(seq[(rval+1):])==0: of.write('F\n')
			else: of.write(line.strip()[(rval+1):]+'\n')
		else:
			of.write(line.strip()+'\n')
	print infile, 'num_partial =', partial

def checkEnd(seq, ada):
	n=len(ada)#last position on adaptor
	head=seq[:n]
	end=ada[n-1]
	rval = -1
	for i in xrange(len(head)):
		if end==head[i] and ada[(n-i-1):n] == head[0:(i+1)]:
			rval=i
	return rval	

def readBlastTab(f):
	trim_index=defaultdict(list)
	for line in f:
		#print line
		parts=line.strip().split()
		try: q, sub, e, qs, qe, ss, se=parts #using fasta
		#try: sub, q, e, ss, se, qs, qe =parts #using blastdb
		except: print 'error', line; sys.exit(1)
		sub, q, ss, se, qs, qe = q, sub, qs, qe, ss, se #use subject db as trim
		try: e = float(e)
		except: e=1
		#if e>10: continue
		#print q
		q=int(q.rsplit('_', 1)[1])
		#print q
		qs, qe=int(qs)-1, int(qe)-1 #change to 0 indexed
		if qs>qe: qs,qe=qe,qs
		trim_index[q].append((qs, qe))
	#print 'len(trim_index)', len(trim_index)
	return trim_index

trim_index = readBlastTab(f2)
f2.close()

i=0 #this is to be consistent with fq2faID.py
num_adaptors, left, right= 0,0,0
for line in f:
	i+=1
	if i%4==1:
		lineno=i/4
		seqid='@s'+str(lineno)+'_'+pair_end+'_'+label
		print >>of, seqid
	elif i%2==0:
		seq=line.strip()
		n=len(seq)
		hits=trim_index[lineno]
		#hits=sorted(hits, key=itemgetter(1), reverse=True)
		minright=n
		maxleft=0
		for hit in hits:
			x,y=hit
			# if i%4==2:
				# print x, y, revcomp(seq[x:y+1])
			mid=(x+y)/2
			if mid>0.5*n: #right
				if x< minright:
					minright=x#seq=seq[0:y]
			else:  #left adaptors
				if y>maxleft:
					maxleft=y
			num_adaptors+=1
		seq=seq[(maxleft+1):minright]
		if len(seq)==0:
			seq=line[0]
		if minright < n: right+=1
		if maxleft  > 0: left+=1
		#if len(hits)>1: print hits, line.strip(); print seq
		print >> of, seq
	else:
		print >>of, line.strip()
print label, '3prime_adaptors = ', right/2 #counted on seq and qseq, so need halfing
print label, '5prime_adaptors = ', left/2 #counted on seq and qseq, so need halfing

f.close()
of.close()

infile = sys.argv[3]+'.tmp' # the file then input to trimend
outfile = sys.argv[3] # output file 
trimAdaptorEnd(infile, outfile, 'CCTTGAAGGCGGACTGTGAG')
