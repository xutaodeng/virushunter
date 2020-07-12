#!/usr/bin/env python

from collections import defaultdict
import operator
import sys
from Bio import pairwise2
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

t = "GTTTCCCACTGGAGGATA"
t1 = revcomp('ACACTCTTTCCCTACACGACGCTCTTCCGATCT')
t2 = revcomp('GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT')
l1 = len(t1) #length of adaptor
l2 = len(t2) #length of adaptor

def trimfqTailB(seqFile, f, of):
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

def trimfqTailS(seqFile, f, of, illumina, phred=13): #illumina 33 or 64
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
			if read.startswith('GGGGG'):
				for r in xrange(len(read)):
					if read[r]!='G':
						s=r+1
						Gcount+=1
						break
			elif read.startswith('CCCCC'):
				for r in xrange(len(read)):
					if read[r]!='C':
						s=r+1
						Ccount+=1
						break
			if read.endswith('GGGGG'):
				for r in xrange(len(read)-1,-1,-1):
					if read[r]!='G':
						t=r
						Gcount+=1
						break
			elif read.endswith('CCCCC'):
				for r in xrange(len(read)-1,-1,-1):
					if read[r]!='C':
						t=r
						Ccount+=1
						break
		elif i%4==3:
			qid=line.strip()
		elif i%4==0: #quality string
			qseq=line.strip()
			index=0
			for x in qseq:
				if ord(x)-illumina <= phred:
					break
				index+=1
			if t<index: index=t
			if index!=0:
				x=len(qseq)
				print >>of, '\n'.join([id, read[s:index],qid,qseq[s:index]])
				counter+=1
				rlen+=(len(qseq)-index)
			else:
				print >>of, '\n'.join([id, read, qid, qseq])
	print seqFile, 'polyC = ', Ccount
	print seqFile, 'polyG = ', Gcount
	#print '\n average tail removed:', float(rlen)/counter, 'nline', i

#return the length of sequence to be trimmed off on read
def len_adaptor(seq, adaptor, l):
	#global t, l
	#print 'seq', seq
	if seq.startswith(t):
		return l
	else:
		try: seq_a, adaptor_a, score, start, end = pairwise2.align.localms(seq, adaptor, 1, 0, -5, -1, one_alignment_only=True)[0]
		except: return 0 #no alignment
	# return the adaptor len on seq
	s_start = start - seq_a.count('-', 0, start)
	a_end = end - adaptor_a.count('-', 0, end)
	#print s_start, a_end, score
	#print seq_a
	#print adaptor_a
	if score>= l-2 or (s_start==0 and a_end==l and score/(end-start)>0.8): #90% matches in alignment
		s_end = end - seq_a.count('-', 0, end)
		#print '****************************************************', s_end
		return s_end
	else:
		return 0

# def trimLen(f, ncpu):
	# i = 0
	# rval=[]
	# for line in f:
		# i+=1
		# if i%4==2: #sequence string
			# rval.append(0)
	# return rval

#return the length of each adaptor in a list
def trimLen2(f, adaptor, length, ncpu):
	i = 0
	pool = Pool(processes=ncpu)
	results =[]
	rval=[]
	bufsize = 10000000
	for line in f:
		i+=1
		seq = line.strip()
		if i%bufsize==0:
			print i
			for result in results:
				rval.append(result.get())
			results=[]
		if i%4==2: #sequence string
			#arg = seq[0:l], adaptor
			results.append(pool.apply_async(len_adaptor, (seq[0:length], adaptor, length)))
	for result in results:
		rval.append(result.get())
	pool.close()
	pool.join()

	return rval


def trimfq(f, of, alens, blens):
	i = 0
	counter=0 #reads counter
	ftrim = 0
	rtrim = 0
	#tt = revcomp("GTTTCCCACTGGAGGATA")
	for line in f:
		i+=1
		seq = line.strip()
		#if i%4000000==0:
		#	print i
		if i%4==2: #sequence string
			counter+=1
			rm_flen = alens[counter-1]
			rm_rlen= blens[counter-1]
			n = max(rm_flen, rm_rlen) #length to remove
			if n>0 and n==rm_flen: ftrim+=1
			elif n>0 and n==rm_rlen: rtrim+=1
			print >>of, seq[n:len(seq)]
		elif i%4==0: #quality string
			print >>of, seq[n:len(seq)]
		else: print >>of, line.strip()
	print 'forward trim', ftrim, 'reverse trim', rtrim, 'No. reads', i/4
	print 'forward \%', 100.0*ftrim/(i/4), 'reverse \%', 100.0*rtrim/(i/4)

def trimfqExact(seqFile, f, of):
	ads =['AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT', 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT', \
	    'CAAGCAGAAGACGGCATACGAGAT']
	revads=[revcomp(ads[0]), revcomp(ads[1]), revcomp(ads[2]) ]
	ads.extend(revads)
	i = 0
	count=0 #reads counter
	outlen=0
	ysum=0
	for line in f:
		i+=1
		if i%4==2: #sequence string
			seq = line.strip()
			left, right=10000000000000, -100000000000000
			x,y = 0, len(seq)
			for a in ads:
				ind = seq.find(a)
				if ind==-1: continue
				if ind<left: left=ind
				if ind+len(a) > right: right=ind+len(a)
			if right > 0: #find a match
				count+=1
				if (left+right) > len(seq):
					y = left
				else:
					x=right
				outlen+=y-x
				ysum+=y
				#print right-left, seq[left:right], len(seq), y-x
			print >>of, seq[x:y]
		elif i%4==0: #quality string
			print >>of, seq[x:y]
		else: print >>of, line.strip()
	print seqFile, 'num_trimed = ', count
	if count==0: count+=1
	print seqFile, 'len_after_trimming = ', float(outlen)/count
	print seqFile, 'pos_adaptor = ', float(ysum)/count

if __name__ == '__main__':
	seqFile = sys.argv[1] #'s_combine_2_sequence.txt'
	outfile =sys.argv[2] #'s_combine_trim_2_sequence.txt'
	illumina=sys.argv[3] #33 or 64 or B
	# ncpu=8
	# f = open(seqFile, 'r')
	# alens = trimLen2(f, t1, l1, ncpu)
	# f.close()
	# f = open(seqFile, 'r')
	# blens = trimLen2(f, t2, l2, ncpu)
	# f.close()
	# print 'alens', len(alens)
	# print 'blens', len(blens)
	of = open('tmp.txt', 'w')
	f = open(seqFile, 'r')
	##trimfq(f, of, alens, blens)
	trimfqExact(seqFile, f, of)
	f.close()
	of.close()
	f = open('tmp.txt', 'r')
	of = open(outfile, 'w')
	#print 'trimming tails'
	if illumina=='B': trimfqTailB(seqFile, f, of) #trim by letter 'B'
	else: trimfqTailS(seqFile, f, of, int(illumina)) # trim by phred score, phred 33
	f.close()
	of.close()
	os.remove('tmp.txt')
