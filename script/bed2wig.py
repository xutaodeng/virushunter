#!/usr/bin/env python

import sys
import os.path
import sys
from collections import defaultdict

def readchrosize():
    chrofile='/mnt/cluster/xdeng/script/ChromSizes/mm9.chrom.sizes'
    f=open(chrofile, 'r')
    chrosize={}
    for line in f:
        chro, size = line.strip().split()
        # chro=chro.upper()
        # if chro.startswith('CHRO'):
            # chro=chro.replace('CHRO', '')
        # elif chro.startswith('CHR'):
            # chro=chro.replace('CHR', '')
        # elif chro.startswith('CH'):
            # chro=chro.replace('CH', '')
        size = int(size)
        chrosize[chro]=size
    f.close()
    print 'chro size', chrosize
    return chrofile, chrosize

def loadsam(filename):
	starts=defaultdict(int)
	ends=defaultdict(int)

	if filename[-3:]==".gz":
		import gzip
		f=gzip.open(filename)
	else: f=open(filename, 'r')
	i=0
	for line in f:
		i+=1
		parts=line.strip().split('\t')
		chro, start, end=parts[0:3]
		start, end = int(start), int(end)

		try:
			starts[chro][start]+=1
			ends[chro][end]+=1
		except:
			starts[chro]=defaultdict(int)
			ends[chro]=defaultdict(int)
			starts[chro][start]=1
			ends[chro][end]=1
			
	f.close()
	return starts, ends

def printwig(outputfile, starts, ends, chrosize):
    allkeys={}
    f=open(outputfile, 'w')
    #print >>f, "track type=bedGraph name="+os.path.basename(outputfile)
    for chro in starts.keys():
        allkeys[chro]=list(set(starts[chro].keys())|set(ends[chro].keys()))
        allkeys[chro].sort()

    for chro in allkeys.keys():
        size=chrosize[chro]
        counter=0
        prevkey=0
        allkey=allkeys[chro]
        start=starts[chro]
        end=ends[chro]
        print >>f, 'variableStep  chrom='+chro
        for i in range(len(allkey)):
            key=allkey[i]
            #print chro, prevkey, key, counter
            if prevkey > size: break 
            if counter!=0: print >> f, prevkey, counter
            prevkey=key
            chrolen=key
            if start.has_key(key): counter+=start[key]
            if end.has_key(key): counter-=end[key]
    f.close()


if __name__ == '__main__':
    inputfiles = sys.argv[1] #bed file
    outputfile = inputfiles+'.wig' #wig file
    chrofile, chrosize=readchrosize()
    starts, ends= loadsam(inputfiles)
    printwig(outputfile, starts, ends, chrosize)
    os.system('wigToBigWig '+outputfile+' '+chrofile+' '+outputfile+'.bw')
