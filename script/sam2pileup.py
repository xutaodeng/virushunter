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
    usage = '''usage: %prog input.sam output.pileup
                generate variant .pileup file from .sam file
            '''
    parser = OptionParser(usage=usage)
    parser.add_option("-p", "--percentage_mutation", dest="percentage", type="float", default=0.1, \
      help="percentage threshold of allele to be called [default: %default]")
    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.error("need input.sam output.pileup!")
    return options.percentage, args

def getNextMD(MDlist, n, deletion):
    rval=[]
    tally=0
    #print 'haha MDlist', MDlist
    if deletion: #return the next deletion
        m = MDlist.popleft()
        while m == '0': m = MDlist.popleft() # cigar'3M3I1M4D89M', MD '3A0^AAAA89',
                                          #the zero is not necessary, so skip it
        #m is '^GA' for example
        if len(m)> n+1: #deletion is in two or more parts
            #seperate m into two parts ^G^A
            MDlist.appendleft('^'+m[(n+1):])
            #print 'MDlist', MDlist
            rval = m[0:(n+1)]
        elif len(m)==n+1:
            rval= m
        else: #cigar 2D, m ^C^T, combine into ^CT
            while len(m)<n+1: #50^C0^T40A0^A0^A20T8 50M2D41M2D29M 
                mnext = MDlist.popleft()
                m = m+mnext[1:]
            rval= m
        #print rval
        return rval


    else: #match or mismatch
        #print 'MDlist', MDlist
        while True:
            m = MDlist.popleft()
            try:
                tally+=int(m)
                if tally < n: rval.append(m)
                else: #split the number into two parts
                    part1, part2 = int(m) - (tally-n), tally-n 
                    rval.append(part1)
                    if part2!=0: MDlist.appendleft(part2)
                    return rval
            except:
                tally+=1
                rval.append(m)
                if tally==n: return rval

def decodeCIGARMD(start, cigar, MD, seq):
    rval=[] # key is Position on reference POR, value is POS (position on read sequence
            # and refletter and mutletter pair
    POR=POS=0
    #subcigar=re.sub('[0-9]+[DNHP]','', cigar) #remove DNHP
    N=re.findall(r'[0-9]+', cigar) #numbers
    L=re.findall(r'[DNHPMIS]', cigar) #letters

    MDlist = deque(re.findall(r'([0-9]+|[ACGTN]|\^[ACGTN]+)', MD))
        
    #print 'MDlist', MDlist
    #print 'cigar', N, L
    # re.findall(r'([0-9ACGTN]|\^[ACGTN]+)', '4^CGA3G2AC^GC')
    #    ['4', '^CGA', '3', 'G', '2', 'A', 'C', '^GC']
    #f1,f2=False,False
    for i in range(len(L)):
        l,n = L[i], int(N[i])
        if l=='D' or l=='M':
            if l=='D': 
                m =getNextMD(MDlist, n, True) #returned string and poped MDlist
                #print 'deletion', m
                #m is a string e.g. ^CAG
                if m[0]!='^':
                    print 'wrong format of MD', MD, m, cigar, seq
                refletter='*'
                mutletter = '-'+str(n)+m[1:] #deletion
                rval.append((POR,POS,refletter, mutletter))
                POR+=n
                f1=True
            else: #match or mimatch
                m=getNextMD(MDlist, n, False) #returned string and poped MDlist
                #print 'match', m
                #now m is a list, e.g. ['A', 23, 'C', 'G']
                for ma in m:
                    try:
                        ma = int (ma)
                        POS+=ma
                        POR+=ma
                    except:
                        refletter=ma
                        mutletter=seq[POS]
                        POS+=1
                        POR+=1
                        rval.append((POR,POS,refletter, mutletter))
        elif l=='I':
            f2=True
            refletter='*'
            mutletter='+'+str(n)+seq[POS:POS+n]
            rval.append((POR,POS,refletter, mutletter))
            POS+=n
        elif l=='H' or l=='N':
            POR+=n
        elif l=='P' or l=='S':
            POR+=n
            POS+=n
    return rval

#first scan to get the mutation dictionary
def loadsam1(filename): # first scan to get the mutation positions
    refkeys=[] #sorted keys for mut
    mut={}
    qua={}
    fCount, rCount=defaultdict(int),defaultdict(int)
    indel=defaultdict(list)
    ref=defaultdict(set)
    if filename[-3:]==".gz":
        import gzip
        f=gzip.open(filename)
    else: f=open(filename, 'r')
    for line in f:
        if line[0]=='@': continue       
        parts=line.strip().split('\t')
        (flag, chro, start, mapq, cigar)=parts[1:6]
        #print 'cigar', cigar
        if chro=='*': continue
        try: seq = parts[9]
        except: print line; continue
        qual=parts[10]
        start=int(start)
        strand='0' # 0 is forward, 1 is reverse
        flagbin=bin(int(flag))
        try: strand =flagbin[-5]
        except: pass
        unmapped='0' # '0' is mapped , 1 is unmapped
        try: unmapped =flagbin[-3]
        except: pass
        if unmapped == '1': continue #unmapped query reads
        tags=parts[11:]
        try:
            MD=[tag for tag in tags if tag.startswith('MD:')][0].split(':')[-1]
            rval=decodeCIGARMD(start, cigar, MD, seq)
        except: rval =[]
        for (POR, POS, refletter, mutletter) in rval:
            if strand =='1': mutletter=mutletter.lower()
            if mutletter[0] in ['+', '-']:
                indel[start+POR-1].append(mutletter)
            else: #not indel
                if not mut.has_key(start+POR-1):
                    mut[start+POR-1]=defaultdict(int)
                    qua[start+POR-1]=defaultdict(list)
                mut[start+POR-1][mutletter]+=1
                try: qua[start+POR-1][mutletter].append(qual[POS-1])
                except:
                    print 'line', line
                    #print 'seq', seq, 'len(seq)', len(seq)
                    #print 'mutletter', mutletter, 'POS', POS, 'qual', qual
                    sys.exit()
            if strand =='1': rCount[start+POR-1]+=1
            else: fCount[start+POR-1]+=1
            ref[start+POR-1].add(refletter.upper())
            #if start+POR-1 == 37879855: print 'rval', rval, refletter
    f.close()
    refkeys=sorted(ref.keys())
   #print 'after one round indel', indel
    #print 'ref', ref[37879855]
    #print 'mutkey', mut[37879855]
    return ref, mut, indel, refkeys, rCount, fCount, qua


#second scan to get all the , and . counts
def loadsam2(filename, refkeys): #second scan to get the , and .
    r2Count, f2Count={},{}#defaultdict(),defaultdict()
    if filename[-3:]==".gz":
        import gzip
        f=gzip.open(filename)
    else: f=open(filename, 'r')
    for line in f:
        if line[0]=='@': continue
        parts=line.strip().split('\t')
        (flag, chro, start)=parts[1:4]
        if chro=='*': continue
        start=int(start)
        seq = parts[9]
        end = start +len(seq)-1
        strand='0' # 0 is forward, 1 is reverse
        flagbin=bin(int(flag))
        try: strand =flagbin[-5]
        except: pass
        unmapped='0' # '0' is mapped , 1 is unmapped
        try: unmapped =flagbin[-3]
        except: pass
        if unmapped == '1': continue #unmapped query reads
        k=bisect_left(refkeys, start)-1
        if k<0: k=0
        while k<len(refkeys) and refkeys[k] <=end:
            pos=refkeys[k]
            if pos >=start:
                if strand =='0':
                    try: f2Count[pos]+=1
                    except: f2Count[pos]=1
                else:
                    try: r2Count[pos]+=1
                    except: r2Count[pos]=1
            k+=1
    f.close()
    #print 'f2count', f2Count
    #print 'r2count', r2Count
    return f2Count, r2Count#now mut contain also keys for indel


def summarizeConsensus(outputfile, ref, mut, indel, rCount, fCount, f2Count, r2Count, qua, percentage, chro):
    filehandle=open(outputfile, 'a') #in append mode
    i=0
    for pos in ref.keys():
        for refletter in ref[pos]:
            if refletter=='N': continue
            try:
                covdict=mut[pos]
                quadict=qua[pos]
            except:
                covdict={}
                quadict={}
            try: indel_list=indel[pos]
            except: indel_list=[]
            try: f2=f2Count[pos]
            except: f2=0
            try: f=fCount[pos]
            except: f=0
            try: r2=r2Count[pos]
            except: r2=0
            try: r=rCount[pos]
            except: r=0
            
            nperiod=f2 - f
            ncomma=r2 - r
            if nperiod<0: nperiod=0
            if ncomma<0: ncomma=0
            coverage=nperiod+ncomma
            count=defaultdict(int)
            for (base, freq) in covdict.items():
                coverage+=freq
                if base=='n' or base=='N': continue
                count[base.upper()]+=freq
            coverage+=len(indel_list)
            #coverage = f2+r2
            threshold=coverage*percentage
            if refletter =='*': cnsletter =indel_list[0][0]
            elif threshold > (f + r):
                continue
            else:
                alleles=[]
                if (nperiod+ncomma) > threshold: alleles.append(refletter)
                alleles.extend([key for key in count.keys() if count[key] > threshold])
                cnsletter=IUPAC_r[frozenset(alleles)]
            #if pos == 37879855: print indel_list, refletter, cnsletter, f, r, threshold
            if refletter != cnsletter:
                covlist=[]
                qualist=[]
                for (base, freq) in covdict.items():
                    covlist.extend([base]*freq)
                    qualist.extend(quadict[base])
                covlist.extend(['.']*nperiod)
                covlist.extend([',']*ncomma)
                covlist.extend(indel_list)
                #qualist don't include qualities for '.,+-' because MD does not have it
                if len(qualist)==0: qualist.append('*')
                print >>filehandle, chro, pos, refletter, cnsletter, '_', '_', '_', coverage, ''.join(covlist), ''.join(qualist)
                i+=1
    return i # 'No. of SNP positions piled-up: ', i

#split sam file into temporary chromosome files 
def splitSAMfile (infile):
    chros=set()
    handles = {}
    tmpfiles=[]
    chroList=[]
    f = open(infile, 'r')
    for line in f:
        if line[0]=='@': continue
        parts=line.strip().split('\t')
        (flag, chro, start)=parts[1:4]
        if chro=='*': continue
        if chro not in chros:
            chros.add(chro)
            tmpfile = chro+'_'+os.path.basename(infile)+'_tmp.sam'
            chroList.append(chro)
            tmpfiles.append(tmpfile)
            handle = open (tmpfile, 'w')
            handles[chro]= handle
        print >> handles[chro], line.strip()
    for chro in handles.keys():
        handles[chro].close()
    return (chroList, tmpfiles)

#only consider the quality of mutant bases, i.e., not '.,'
if __name__ == "__main__":
    (percentage, args) = options()
    print 'piling up...'
    infile, outputfile = args
    if os.path.exists (os.path.abspath(outputfile)): 
        os.remove(os.path.abspath(outputfile))

    chroList, inputfiles=splitSAMfile(infile)
    numSNPs=0
    for i in xrange (len(chroList)):
        inputfile = inputfiles[i]
        chro=chroList[i]
        ref, mut, indel, refkeys, rCount, fCount, qua = loadsam1(inputfile)
        f2Count, r2Count = loadsam2(inputfile, refkeys)
        numSNPs += summarizeConsensus(outputfile, ref, mut, indel, rCount, fCount, f2Count, r2Count, qua, percentage, chro)
        os.remove(os.path.abspath(inputfile))
    print "Number of SNPs:", numSNPs
