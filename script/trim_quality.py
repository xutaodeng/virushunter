#!/usr/bin/env python

from collections import defaultdict
import operator
import sys
from multiprocessing import Pool
import os

def revcomp(tag):
	rval = []
	for base in tag:
		if base == 'A': a='T'
		elif base == 'C': a='G'
		elif base == 'G': a='C'
		elif base == 'T': a='A'
		rval.append(a)
	rval.reverse()
	return ''.join(rval)

def trimfqTailB(label, seqFile, f, of):
	i = 0
	rlen=0
	counter=0
	s,t=None, None
	Ccount, Gcount =0,0
	for line in f:
		i+=1
		#if i%1000000==0: print i,' ',
		if i%4==1:
			id=line.strip()
		elif i%4==2: #sequence string
			read=line.strip()
			s,t=0, len(read)
			if read.startswith('GGGG'):
				for r in xrange(len(read)):
					if read[r]!='G':
						s=r+1
						Gcount+=1
						break
			elif read.startswith('CCCC'):
				for r in xrange(len(read)):
					if read[r]!='C':
						s=r+1
						Ccount+=1
						break
			if read.endswith('GGGG'):
				for r in xrange(len(read)-1,-1,-1):
					if read[r]!='G':
						t=r
						Gcount+=1
						break
			elif read.endswith('CCCC'):
				for r in xrange(len(read)-1,-1,-1):
					if read[r]!='C':
						t=r
						Ccount+=1
						break
		elif i%4==3:
			qid=line.strip()
		elif i%4==0: #quality string
			qseq=line.strip()
			try: index = qseq.index('B')
			except: index=len(qseq)
			if t<index: index=t
			x=len(qseq)
			print >>of, '\n'.join([id, read[s:index],qid,qseq[s:index]])
			counter+=1
			rlen+=(len(qseq)-index)
	print seqFile, 'polyC = ', Ccount
	print seqFile, 'polyG = ', Gcount
	#print '\n average tail removed:', float(rlen)/counter, 'nline', i

def trimfqTailS(label, seqFile, f, of, illumina, phred=10): #illumina 33 or 64
	i = 0
	rlen=0
	counter=0
	s,t=None, None
	Ccount, Gcount =0,0
	for line in f:
		i+=1
		#if i%1000000==0: print i,' ',
		if i%4==1:
			id=line.strip()
		elif i%4==2: #sequence string
			read=line.strip()
			s,t=0, len(read)
			#trim 5' polyN
			r=0 #position
			while r< (t-1):
				if read[r]==read[r+1]:
					r+=1
				else:
					break
			if r>0: s=r+1
			# if read.startswith('GGG'):
				# for r in xrange(len(read)):
					# if read[r]!='G':
						# s=r+1
						# Gcount+=1
						# break
			# elif read.startswith('CCC'):
				# for r in xrange(len(read)):
					# if read[r]!='C':
						# s=r+1
						# Ccount+=1
						# break
			# if read.endswith('GGG'):
				# for r in xrange(len(read)-1,-1,-1):
					# if read[r]!='G':
						# t=r
						# Gcount+=1
						# break
			# elif read.endswith('CCC'):
				# for r in xrange(len(read)-1,-1,-1):
					# if read[r]!='C':
						# t=r
						# Ccount+=1
						# break
		elif i%4==3:
			qid=line.strip()
		elif i%4==0: #quality string
			qseq=line.strip()
			index=0
			zz=0
			for x in qseq:
				zz+=1
				if zz<40: pass #do not check front of read
				elif ord(x)-illumina <= phred:
					break
				index+=1
			#if t<index: index=t
			if index!=0 and index>s:
				x=len(qseq)
				print >>of, '\n'.join([id, read[s:index],qid,qseq[s:index]])
				counter+=1
				rlen+=(len(qseq)-index)
			else:
				print >>of, '\n'.join([id, read, qid, qseq])
	# print label, 'polyC = ', Ccount
	# print label, 'polyG = ', Gcount
	print label, 'num_tail_removed = ', counter
	if counter==0: 
		out=0
	else:
		out=float(rlen)/counter
	print label, 'tail_removed_average = ', out

if __name__ == '__main__':
	seqFile = sys.argv[1] #'s_combine_2_sequence.txt'
	outfile =sys.argv[2] #'s_combine_trim_2_sequence.txt'
	illumina=sys.argv[3] #33 or 64 or B
	label=sys.argv[4] #label of input file library
	try: phred= int(sys.argv[5])
	except: phred=10
	f = open(seqFile, 'r')
	of = open(outfile, 'w')
	if illumina=='B': trimfqTailB(label, seqFile, f, of) #trim by letter 'B'
	else: trimfqTailS(label, seqFile, f, of, int(illumina), phred) # trim by phred score, phred 33
	f.close()
	of.close()