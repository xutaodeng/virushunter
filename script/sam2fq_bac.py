#!/usr/bin/env python

import sys
import os
import os.path
from optparse import OptionParser
from collections import defaultdict, deque
from bisect import bisect_left
import re
from operator import itemgetter

# def checkBacSAM(samindex, start, end): # first scan to get the mutation positions
	# fs1, fs2=[],[]
	# for j in range(start, end+1):
		# filename=samindex+'_1'+'_'+str(j)+'.sam'
		# fs1.append(open(filename, 'r'))
	# for j in range(start, end+1):
		# filename=samindex+'_2'+'_'+str(j)+'.sam'
		# fs2.append(open(filename, 'r'))

	# fs = fs1+fs2
	# for f in fs:
		# count+=1
		# line1 = f1.readline() 
		# line2 = f2.readline() 
		# if not line1: end=True; break
	# if end: break
	# for f1,f2 in zip(fs1, fs2):
		# f1.close()
		# f2.close()
	# print filename, 'percentHuman = ', 100*(1-(float(bacCount)/count))

def processBacSAM(samindex, start, end, outfile,outfile2): # first scan to get the mutation positions
	fs1, fs2=[],[]
	for j in range(start, end+1):
		filename=samindex+'_1'+'_'+str(j)+'.sam'
		fs1.append(open(filename, 'r'))
	for j in range(start, end+1):
		filename=samindex+'_2'+'_'+str(j)+'.sam'
		fs2.append(open(filename, 'r'))
	of=open(outfile, 'w')
	of2=open(outfile2, 'w')
	count=0
	bacCount,bacCount1,bacCount2=0,0,0
	humanCount,humanCount1,humanCount2=0,0,0
	end = False
	while 1:
		lines =[]
		bac1,bac2 = False,False
		human1,human2 = False,False
		index =0
		for f1,f2 in zip(fs1, fs2):
			line1 = f1.readline() 
			line2 = f2.readline() 
			if not line1: end=True; break
			parts=line1.strip().split('\t')
			chro, seq=parts[2], parts[9]
			if chro !='*' and len(seq)>20:
				if start==0 and index ==0:
					human1=True
				else:
					bac1=True
			parts=line2.strip().split('\t')
			chro, seq=parts[2], parts[9]
			if chro !='*' and len(seq)>20:
				if start==0 and index ==0:
					human2=True
				else:
					bac2=True
			#if line1==line2: print 'haha'
			index+=1
		if end: break
		count+=1
		if bac1: bacCount1+=1
		if bac2: bacCount2+=1
		if bac1 and bac2:  bacCount+=1; continue
		if human1: humanCount1+=1
		if human2: humanCount2+=1
		if human1 and human2:  humanCount+=1; continue
		parts=line1.strip().split('\t')
		parts2=line2.strip().split('\t')
		(name, flag, chro, sstart, mapq, cigar)=parts[0:6]
		(name2, flag2, chro2, start2, mapq2, cigar2)=parts2[0:6]
		name='@'+name
		name2='@'+name2
		try: seq = parts[9]
		except: print line; seq='A'
		try: seq2 = parts2[9]
		except: print line2; seq2='A'
		qual=parts[10]
		qual2=parts2[10]
		if bac1: seq='A'; qual='G'
		if bac2: seq2='A'; qual2='G'
		print >>of, '\n'.join([name, seq, '+',qual])
		print >>of2, '\n'.join([name2, seq2, '+',qual2])

	of.close()
	of2.close()
	for f1,f2 in zip(fs1, fs2):
		f1.close()
		f2.close()

	print filename, 'percentHuman = ', 100*float(humanCount)/count
	print filename, 'percentBac = ', 100*float(bacCount)/count #, bacCount1, bacCount2, bacCount, count
	#print filename, 'filtered', filter, 'unmapped', un, 'DNA', geno, 'mRNA', ma, \
	#	  'HumanPercent', float(ma+geno)/(un+ma+geno)*100, '%'

def processSingleBacSAM(samindex, start, end, outfile): # first scan to get the mutation positions
	fs1=[]
	for j in range(start, end+1):
		filename=samindex+'_1'+'_'+str(j)+'.sam'
		fs1.append(open(filename, 'r'))
	of=open(outfile, 'w')
	count=0
	bacCount, humanCount=0,0
	end = False
	while 1:
		lines =[]
		bac=False
		human=False
		index =0
		for f1 in fs1:
			line1 = f1.readline() 
			if not line1: end=True; break
			parts=line1.strip().split('\t')
			chro, seq=parts[2], parts[9]
			if chro !='*' and len(seq)>20:
				if start==0 and index ==0:
					human=True
				else:
					bac=True
			index+=1
		count+=1
		if end: break
		if bac: bacCount+=1; continue
		if human: humanCount+=1; continue
		parts=line1.strip().split('\t')
		(name, flag, chro, sstart, mapq, cigar)=parts[0:6]
		name='@'+name
		try: seq = parts[9]
		except: print line; seq='A'
		qual=parts[10]
		print >>of, '\n'.join([name, seq, '+',qual])
	for f1 in fs1:
		f1.close()
	of.close()
	print filename, 'percentBac = ', 100*float(bacCount)/count
	print filename, 'percentHuman = ', 100*float(humanCount)/count


if __name__ == "__main__":
	try: processBacSAM(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5]) #pair end
	except: processSingleBacSAM(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]) #single end
