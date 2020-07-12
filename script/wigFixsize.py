#!/usr/bin/env python

import sys
import os.path
import re
from optparse import OptionParser
from collections import defaultdict

def options():
    usage = '''usage: %prog organism input.wig output.wig output.bw
                fix input.wig to chromo sizes and generate bigwig
                example: wigFixsize.py hg18 input.wig output.wig output.bw
                support hg18, hg19, mm10
            '''

def readchrosize(organism):
    directory=os.path.dirname(sys.argv[0])
    chrofile=directory+'/ChromSizes/'+organism+'.chrom.sizes'
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
    return chrosize

def printwig(inwig, outwig, outbw, chrosize, organism):
    f=open(inwig, 'r')
    of=open(outwig, 'w')
    print 'fixing wig'
    #of.write(f.readline().strip()+'\n')#header
    for line in f:
        parts=line.strip().split()
        if len(parts)==3:
            chro=parts[1].split('=')[0]
            print >> of, line.strip()
            print line
            continue
        start=''
        if len(parts)==2: start, counter = parts
        else: continue
        try: start= int(start)
        except: print line; sys.exit()
        try: size=chrosize[chro]
        except: 
            size=10000000000000
            #print 'warning: did ot find chrosize for', ch, 'using', size
            if start >=size: continue
            print >> of, start, counter
    f.close()
    of.close()
#wigToBigWig in.wig chrom.sizes out.bw
    directory=os.path.dirname(sys.argv[0])
    chrofile= directory+'/ChromSizes/'+organism+'.chrom.sizes'
    print 'wigToBigwig'
    os.system('wigToBigWig '+outwig+' '+chrofile+' '+outbw)

if __name__ == '__main__':
    options()
    try: organism , inwig, outwig, outbw= sys.argv[1:5]
    except: 
        organism , inwig = sys.argv[1:3]
        outwig, outbw= inwig+'2', inwig+'.bw'
    chrosize = readchrosize(organism)
    printwig(inwig, outwig, outbw, chrosize, organism)
