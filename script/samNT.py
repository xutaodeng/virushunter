#!/usr/bin/env python

import sys
import os
import os.path
from optparse import OptionParser
from collections import defaultdict, deque
from bisect import bisect_left
import re
from operator import itemgetter
import linecache
import re
import os.path



def CacheLines(fname): 
	cache={}
	f = open(fname, 'r')
	i=0
	start, end = 0,0
	header= None
	for line in f:
		i+=1
		if line.strip().startswith('>'):
			end=i-1
			if header!=None: cache[header] = (start, end)
			header = line.strip()[1:]
			start=i+1
	if header!=None: cache[header] = (start, i)
	return cache

def getSeq(cachename, cache, header):
	seq=[]
	start, end = cache[header]
	for i in xrange(start, end+1):
		seq.append(linecache.getline(cachename, i).strip())
	return ''.join(seq)

def processSAM(key, wd, base, startInd, endInd, cahche): # first scan to get the mutation positions
	countfile=wd+'/'+base+'/clark/'+key+'.count'
	of=open(countfile,'w')
	of2=open(countfile+'.csv','w')

	fs1=[]
	for j in range(startInd, endInd+1):
		filename=wd+'/NT/'+ key+'_'+str(j)+'.sam'
		fs1.append(open(filename, 'r'))

	count1=defaultdict(int)
	counts2=defaultdict(list)

	end = False
	while 1:
		hit=False
		for f1 in fs1:
			line1 = f1.readline() 
			if not line1: end=True; break
			parts=line1.strip().split('\t')
			name, chro, seq=parts[0], parts[2], parts[9]
			if chro !='*' and len(seq)>20:
				chro=chro.replace(',','')
				cat, clas, fam, species=chro.split('$')
				#if ',' in species: print chro
				count1[cat]+=1
				category2=chro.replace('$', ',')
				counts2[category2].append(name)
				hit=True
				break
		if not hit:
			#cat, clas, fam, species=chro.split('$')
			count1['NA']+=1
			#count2['NA$NA$NA$NA']+=1

		if end: break

		# parts=line1.strip().split('\t')
		# (name, flag, chro, sstart, mapq, cigar)=parts[0:6]
		# name='@'+name
		# try: seq = parts[9]
		# except: print line; seq='A'
		# qual=parts[10]
		# print >>of, '\n'.join([name, seq, '+',qual])
	for f1 in fs1:
		f1.close()

	for key in count1.keys():
		of.write(key+'\t'+str(int(count1[key]))+'\n')
	of.close()

	of2.write('category,class,family,species,count\n')

	x=[(key, len(snames)) for (key, snames) in counts2.items()]
	sorted_x = sorted(x, key=itemgetter(1), reverse=True)

	for key, val in sorted_x:
		seqnames=counts2[key]
		of2.write(key+','+str(int(val))+'\n')
		key2=key.replace(',', '_').replace('/', '_')
		outfa=wd+'/'+base+'/clark/fasta/'+os.path.basename(countfile)+'.csv.'+key2+'.fa'
		fa=open(outfa, 'w')
		for seqname in seqnames:
			query_nt = getSeq(fafile, cache, seqname)
			fa.write('>'+seqname+'\n')
			fa.write(query_nt+'\n')
		fa.close()

	of2.close()

if __name__ == "__main__":
	key=sys.argv[1]
	wd=sys.argv[2]
	base=sys.argv[3]
	startInd=int(sys.argv[4])
	endInd=int(sys.argv[5])
	
	fafile=wd+'/fastq/'+key+'.fa'
	cache = CacheLines(fafile)
	processSAM(key, wd, base, startInd, endInd, cache) #single end
