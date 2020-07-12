#!/usr/bin/env python

import sys
import os
import os.path
from optparse import OptionParser
from collections import defaultdict, deque
from bisect import bisect_left
from codon import *
import re
from operator import itemgetter

def options():
    usage = '''usage: %prog input.sam output.sam
                remove unmatched reads to minimize file size
            '''
    parser = OptionParser(usage=usage)
    parser.add_option("-r", "--header_remove",
                  action="store_true", dest="hr", default=False,
                  help="remove header in the output")
    parser.add_option("-u", "--unalign_remove",
                  action="store_true", dest="ur", default=False,
                  help="remove unaligned records in the output")
    parser.add_option("-g", "--gzip",
                  action="store_true", dest="gzip", default=False,
                  help="gzip output")   
    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.error("need input.sam output.sam!")
    return options.hr, options.ur, options.gzip, args


#first scan to get the mutation dictionary
def processsam(infile, outfile, hr, ur, gz): # first scan to get the mutation positions
    if infile[-3:]==".gz":
        import gzip
        f=gzip.open(infile)    
    else: f=open(infile, 'r')

    if gz==True:
        import gzip
        of=gzip.open(outfile, 'wb')    
    else: of=open(outfile, 'w')

    i=0
    for line in f:
        i+=1
        if i%100000==0: print 'line', i
        if line[0]=='@' and hr == True: continue
        if line[0]!='@':
            parts=line.strip().split('\t')
            (flag, chro, start, mapq, cigar)=parts[1:6]
            if ur == True and chro=='*': continue
        of.write(line)

    f.close()
    of.close()

#only consider the quality of mutant bases, i.e., not '.,'
if __name__ == "__main__":
    hr, ur, gz, args = options()
    infile, outfile = args
    processsam(infile, outfile, hr, ur, gz)

