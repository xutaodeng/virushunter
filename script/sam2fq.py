#!/usr/bin/env python

import sys
import os
import os.path
from optparse import OptionParser
from collections import defaultdict, deque
from bisect import bisect_left
import re
from operator import itemgetter

genomic=set(['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'MT'])
virus=set(['unclass', 'Viruses']) #unclassified NT sequences
want=set(['unclass', 'Viruses', 'unmap'])

def CIGARMatch(cigar):
#first scan to get the mutation dictionary
	N=re.findall(r'[0-9]+', cigar) #numbers
	L=re.findall(r'[DNHPMIS]', cigar) #letters
	for i in range(len(L)):
		l,n = L[i], int(N[i])
		if l=='M' and int(n)>=50:
			return True
	return False
def getType(lines):
	align= 'comment'
	oldname, name, seq, qual=None, None, None, None
	i=0
	for line in lines:
		i+=1
		#print i, line
		if line[0]=='@': break
		parts=line.strip().split('\t')
		(name, flag, chro, start, mapq, cigar)=parts[0:6]
		name='@'+name
		if oldname!=None and name!=oldname: 
			print 'oldname', oldname, 'name', name
			for ln in lines:
				print ln.strip().split()[0]
			sys.exit(1)
		oldname=name
		try: seq = parts[9]
		except: seq='A'
		qual=parts[10]
		start=int(start)
		if len(seq)==1:
			align='filter'
			break
			#print >>of, '\n'.join([name, 'A', '+','A'])
		elif chro == '*':
			align='unmap' #keep looking
		elif chro =='Viruses':
			align=chro
			break
		elif CIGARMatch(cigar): #match non virus seq in NT
			align=chro
			break
		else: #match non virus but not long enough
			align='unmap' #this is POSSIBLY Virus, keep looking
	#print align, name, seq, qual
	#print len(lines)
	return align, name, seq, qual

# def processNTSAMs(samindex, outfile): # first scan to get the mutation positions
	# pie=defaultdict(int)
	# fs=[]
	# for j in range(1, 13):
		# filename=samindex+'_'+str(j)+'.sam'
		# fs.append(open(filename, 'r'))
	# of=open(outfile, 'w')
	# end=False
	# k=0
	# while 1:
		# lines=[]
		# for f in fs:
			# line = f.readline() 
			# if not line: end=True; break
			# lines.append(line)
		# if end: break
		# align, name, seq, qual=getType(lines)
		# k+=1
		# #if k%10000==0: print 'k',k, name
		
		# pie[align]+=1
		# if align=='comment': continue
		# if align in want:
			# print >>of, '\n'.join([name, seq, '+',qual])
		# else:
			# print >>of, '\n'.join([name, 'A', '+','A'])
	# of.close()
	# for f in fs:
		# f.close()
	# sum=0
	# for key in pie.keys():
		# print key, pie[key]
		# if key!='filter': continue
			# sum+=pie[key]
	# for key in pie:
		# print samindex, key, '=', pie[key]
	# print samindex, 'percent_candidate_virus =', 100.0*(pie['unclass']+ pie['Viruses']+pie['unmap'])/sum

	
# def processNTSAM(filename, outfile): # first scan to get the mutation positions
	# filter=0
	# pie=defaultdict(int)
	# if filename[-3:]==".gz":
		# import gzip
		# f=gzip.open(filename)
	# else: f=open(filename, 'r')
	# of=open(outfile, 'w')
	# for line in f:
		# if line[0]=='@': continue	   
		# parts=line.strip().split('\t')
		# (name, flag, chro, start, mapq, cigar)=parts[0:6]
		# name='@'+name
		# try: seq = parts[9]
		# except: seq='A'
		# qual=parts[10]
		# start=int(start)
		# if len(seq)==1:
			# filter+=1
			# print >>of, '\n'.join([name, 'A', '+','A'])
		# elif chro == '*':
			# print >>of, '\n'.join([name, seq, '+',qual])
			# pie[chro]+=1
		# elif chro in virus:
			# print >>of, '\n'.join([name, seq, '+',qual])
			# pie[chro]+=1
		# else:
			# if CIGARMatch(cigar): #match non virus seq in NT
				# print >>of, '\n'.join([name, 'A', '+','A'])
				# pie[chro]+=1
			# else: #match non virus but not long enough
				# print >>of, '\n'.join([name, seq, '+',qual])
				# pie['*']+=1
	# f.close()
	# print filename, 'filtered', filter
	# sum=0
	# for key in pie.keys():
		# print key, pie[key]
		# sum+=pie[key]
	# print 'candidate virus percent', 100.0*(pie['unclass']+ pie['Viruses'])/sum+'%'
	
def processHumanSAM(filename, filename2, outfile,outfile2): # first scan to get the mutation positions
	un, ma, geno, filter=0,0,0,0
	if filename[-3:]==".gz":
		import gzip
		f=gzip.open(filename)
	else: f=open(filename, 'r')
	of=open(outfile, 'w')
	if filename2[-3:]==".gz":
		import gzip
		f2=gzip.open(filename2)
	else: f2=open(filename2, 'r')
	of2=open(outfile2, 'w')
	count=0
	for line in f:
		line2=f2.readline()
		if line[0]=='@': continue
		count+=1
		parts=line.strip().split('\t')
		parts2=line2.strip().split('\t')
		(name, flag, chro, start, mapq, cigar)=parts[0:6]
		(name2, flag2, chro2, start2, mapq2, cigar2)=parts2[0:6]
		name='@'+name
		name2='@'+name2
		#print 'cigar', cigar
		try: seq = parts[9]
		except: print line; seq='A'
		try: seq2 = parts2[9]
		except: print line2; seq2='A'
		qual=parts[10]
		qual2=parts2[10]
		# strand='0' # 0 is forward, 1 is reverse
		# flagbin=bin(int(flag))
		# unmapped='0' # '0' is mapped , 1 is unmapped
		# try: unmapped =flagbin[-3]
		# except: pass
		# if len(seq)==1:
			# filter+=1
			# print >>of, '\n'.join([name, 'A', '+','A'])
		if chro == '*' and len(seq)>20 and chro2 == '*' and len(seq2)>20: 
			un+=1
			print >>of, '\n'.join([name, seq, '+',qual])
			print >>of2, '\n'.join([name2, seq2, '+',qual2])
		# elif chro in genomic:
			# geno+=1
			# print >>of, '\n'.join([name, 'A', '+','A'])
		# else:
			# ma+=1
			# print >>of, '\n'.join([name, 'A', '+','A'])
	f.close()
	print filename, 'percentHuman = ', 100*(1-(float(un)/count))
	#print filename, 'filtered', filter, 'unmapped', un, 'DNA', geno, 'mRNA', ma, \
	#	  'HumanPercent', float(ma+geno)/(un+ma+geno)*100, '%'

def processHumanSAM_single(filename, outfile): # first scan to get the mutation positions
	un, ma, geno, filter=0,0,0,0
	if filename[-3:]==".gz":
		import gzip
		f=gzip.open(filename)
	else: f=open(filename, 'r')
	of=open(outfile, 'w')

	count=0
	for line in f:
		if line[0]=='@': continue
		count+=1
		parts=line.strip().split('\t')
		(name, flag, chro, start, mapq, cigar)=parts[0:6]
		name='@'+name
		#print 'cigar', cigar
		try: seq = parts[9]
		except: print line; seq='A'
		qual=parts[10]
		if chro == '*' and len(seq)>20: 
			un+=1
			print >>of, '\n'.join([name, seq, '+',qual])
	f.close()
	print filename, 'percentHuman = ', 100*(1-(float(un)/count))
	
#only consider the quality of mutant bases, i.e., not '.,'
if __name__ == "__main__":
	#samindex, outfile =sys.argv[1], sys.argv[2]
	#processNTSAM(filename, outfile)
	#processNTSAMs(samindex, outfile)
	try: processHumanSAM(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]) #pair end
	except: processHumanSAM_single(sys.argv[1], sys.argv[2]) #single end
	#processHumanSAM(filename, outfile)
