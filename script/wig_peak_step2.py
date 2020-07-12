#!/usr/bin/env python

import os
import os.path
import sys
from collections import defaultdict
from operator import itemgetter
import bisect

rval2={}
rval=defaultdict(float)

# def get_wig(peaks, wig):
	# global rval
	# for start, end in peaks:
		# center=(start+end)/2
		# for wstart, value in wig:
			# wcenter = wstart+span/2
			# if start < wcenter < end:
				# dist= (center-wcenter)/span
				# rval[dist]+=value
def get_wig(chro, peaks, wig):
	global rval, rval2
	#wig=sorted(wig,key=itemgetter(0)) #sort by starts
	#ws = [x for (x, y) in wig] #wig starts sorted
	i=0
	j=0
	if len(peaks)==0: return
	
	while 1:
		start, end  = peaks[i]
		peakname= chro+'_'+str(start)+'_'+str(end)
		if not rval2.has_key(peakname):
			rval2[peakname]=defaultdict(float)
		center=(start+end)/2
		#j=bisect.bisect_left(ws, start)-1
		#if j<0: j=0
		wstart, value = wig[j]
		wcenter = wstart+span/2
		if start < wcenter < end:
			dist= (center-wcenter)/span
			rval2[peakname][dist]+=value
			rval[dist]+=value
		elif wcenter >= end:
			i+=1
		j+=1
		if i >=len(peaks) or j>=len(wig): break

span=10
f1=open(sys.argv[1], 'r')
f2=open(sys.argv[2], 'r')
of=open(sys.argv[3], 'w')
of2=open(sys.argv[3]+'.matrix', 'w')
of3=open(sys.argv[3]+'.log', 'w')

print 'loading peak'
peaks=defaultdict(list)
for line in f2:
	parts=line.strip().split()
	chro, start, end = parts[0:3]
	start, end = int(start), int(end)
	peakcenter=(start+end)/2
	peaks[chro].append((peakcenter-2000,peakcenter+2000))
for chro in peaks.keys():
	peaks[chro].sort(key=itemgetter(0))
	#print chro, 'num peaks', len(peaks[chro])

wig=[]
for line in f1:
	parts=line.strip().split()
	if len(parts)==3: 
		if wig!=[]: 
			print 'getting peaks'
			get_wig(chro, peaks[chro], wig)
		wig=[]
		chro=parts[1].split('=')[1]
		print 'loading wig', chro
		of3.write('loading wig '+chro+'\n')
		span=int(parts[2].split('=')[1])
	elif len(parts)>3: continue
	else:
		start, value = parts
		wig.append((int(start), float(value)))

get_wig(chro, peaks[chro], wig)

print 'len(wig)', len(wig)
print 'rval.keys()', rval.keys()
ks=sorted(rval.keys())
for k in ks:
	of.write(str(k)+' '+str(rval[k])+'\n')

of2.write('UNIQID')
for k in ks:
	of2.write('\t'+str(k))
of2.write('\n')

for pk in rval2.keys(): #peaknames
	i=0
	for k in ks: 
		if i==0:
			of2.write(pk)
		of2.write('\t'+str(rval2[pk][k]))
		i+=1
	of2.write('\n')

of.close()
of2.close()
of3.close()
