#!/usr/bin/env python
from collections import defaultdict
import operator
import sys
import os.path
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
    
def openBarcodes(outdir, ff, pair):
    tags={}
    fhandles={}
    for line in ff:
        parts = line.strip().split()
        id, tag = parts[0:-1], parts[-1]
        id = '_'.join(id)
        id = id.replace('/', '_')
        tags[id]= (tag, revcomp(tag))
        if pair: #paried files
            id1 = id + '_1'
            id2 = id + '_2'
            fhandles[id1]=open(os.path.join(outdir, id1+'_sequence.txt'), 'w')
            fhandles[id2]=open(os.path.join(outdir,id2+'_sequence.txt'), 'w')
        else: #single end file
            fhandles[id]=open(os.path.join(outdir, id+'_sequence.txt'), 'w')
    if pair:
        fhandles['no_tag_1']=open(os.path.join(outdir, 'no_tag_1_sequence.txt'), 'w')
        fhandles['no_tag_2']=open(os.path.join(outdir, 'no_tag_2_sequence.txt'), 'w')
    else:
        fhandles['no_tag']=open(os.path.join(outdir, 'no_tag_sequence.txt'), 'w')
    return tags, fhandles

def closeBarcodes(ff, fhandles):
    for (id, handle) in fhandles.items():
        handle.close()
    ff.close()

#process single fastq file
def processfq1(f):
    global fhandles, miss, hit
    i = 0
    for line in f:
        i+=1
        if i%4==1:
            header=line.strip()
            id = processHeader(header)
            if id =='no_tag': miss+=1
            else: hit+=1
        elif i%4==2:
            seq = line.strip()
        elif i%4==3:
            quaID = line.strip()
        elif i%4==0 and i>1:
            qual = line.strip()
            fhandles[id].write(header+'\n'+seq+'\n'+quaID+'\n'+qual+'\n')

#process pair-end fastq file
#use the barcode of file 1 to debarcode
def processfq2(f1, f2):
    global fhandles, miss, hit
    i = 0
    line1 = 'dummy'
    while line1 != "":
        line1 = f1.readline()
        line2 = f2.readline()
        i+=1
        if i%4000000==0: print i/4
        if i%4==1:
            header1=line1.strip()
            header2=line2.strip()
            id = processHeader(header1)
            if id =='no_tag': miss+=1
            else: hit+=1
            id1 = id +'_1'
            id2 = id +'_2'
        elif i%4==2:
            seq1 = line1.strip()
            seq2 = line2.strip()
        elif i%4==3:
            quaID1 = line1.strip()
            quaID2 = line2.strip()
        elif i%4==0 and i>1:
            qual1 = line1.strip()
            qual2 = line2.strip()
            fhandles[id1].write(header1+'\n'+seq1+'\n'+quaID1+'\n'+qual1+'\n')
            fhandles[id2].write(header2+'\n'+seq2+'\n'+quaID2+'\n'+qual2+'\n')
   
def processHeader(header):
    global tags
    for (id, tag) in tags.items():
        ftag, rtag = tag[0], tag[1]
        #if ftag in header: 
        #    return id
        if rtag in header:
            return id
    return 'no_tag'
        

seqFile1 = sys.argv[1] #'s_8_1_trim_sequence.txt'
seqFile2 = sys.argv[2] #'s_8_2_trim_sequence.txt'
pair = True
outdir = 'fastq'
if not os.path.exists(outdir): os.mkdir(outdir)
barcodeFile = sys.argv[3] #'ILL_barcode.txt'
f1 = open(seqFile1, 'r')
f2 = open(seqFile2, 'r')
ff = open(barcodeFile, 'r')

tags, fhandles = openBarcodes(outdir, ff, pair)
print tags

hit, miss = 0, 0
processfq2(f1, f2)

closeBarcodes(ff, fhandles)
f1.close()
f2.close()
print 'hit', hit, 'miss', miss, 'hit_rate', float(hit)/(miss+hit)
