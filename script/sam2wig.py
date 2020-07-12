#!/usr/bin/env python

import sys
import os.path
import re
from optparse import OptionParser
from collections import defaultdict
import uuid


def options():
    usage = '''usage: %prog -e 0 -n 1000000 input1.sam input2.sam ... output.wig
                convert one .sam file to .wig file with extention of xxbp added 
                or aggregate multiple .sam file to .wig file
            '''
        
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--strand", dest="strand", default="both", \
      help="strand: forward, reverse, both [default: %default]")
    # parser.add_option("-o", "--organism", dest="organism", default="hg19", \
      # help="choose hg19, hg18, or mm10: [default: %default]")
    parser.add_option("-n", "--normalize", dest="normalize", default="0", \
      help="choose normalize to a fixed number of reads such as 1000000: [default: %default]")
    parser.add_option("-e", "--extension", dest="extension", default=0, \
      help="extension: extend the read length by this much [default: %default]")
    (options, args) = parser.parse_args()

    return (options.strand, options.extension,  options.normalize, args)

# def readchrosize(organism):
    # directory=os.path.dirname(sys.argv[0])
    # chrofile=directory+'/ChromSizes/'+organism+'.chrom.sizes'
    # # print 'chrosize file', chrofile
    # # f=open(chrofile, 'r')
    # # chrosize={}
    # # for line in f:
        # # chro, size = line.strip().split()
        # # chro=chro.upper()
        # # if chro.startswith('CHRO'):
            # # chro=chro.replace('CHRO', '')
        # # elif chro.startswith('CHR'):
            # # chro=chro.replace('CHR', '')
        # # elif chro.startswith('CH'):
            # # chro=chro.replace('CH', '')
        # # size = int(size)
        # # chrosize[chro]=size
    # # f.close()
    # # print 'chro size', chrosize
    # return chrosize, chrofile

def loadsam(filenames, strand, extension=0):
    starts=defaultdict(int)
    ends=defaultdict(int)
    unaligned, aligned, spliced, forward, reverse = 0,0,0,0,0 #bad and good alignment
    #global unaligned
    #global aligned
    gsize, bases=0, 0
    for filename in filenames:
        if filename[-3:]==".gz":
            import gzip
            f=gzip.open(filename)
        else: f=open(filename, 'r')
        i=0
        for line in f:
            i+=1
            if i%100000==0: print >> sys.stderr, i, 'aligned', aligned
            if line[0]=='@' or line[0]=='#':
                if line.startswith('@SQ'):
                    gsize+=int(line.strip().split()[-1].split(':')[-1])
                continue
            parts=line.strip().split('\t')
            try: (flag, chro, start)=parts[1:4]
            except:
                #print 'sam format error, ignored line' +str(i)+': '+line
                continue
            try: flagbin=bin(int(flag))
            except: continue
            unmapped='0' # '0' is mapped , 1 is unmapped
            try: unmapped =flagbin[-3]
            except: pass
            if unmapped == '1' or chro=='*':
                unaligned+=1
                continue
            aligned+=1
            try: stran =flagbin[-5]
            except: stran = '0'
            if stran =='1': reverse+=1
            else: forward+=1
            if strand == 'forward' and stran != '0': continue
            elif strand == 'reverse' and stran !='1': continue
            #print 'stran', stran, flagbin
            cigar = parts[5]
            
            start= long(start)
            N=[]
            if cigar.find('N') == -1: # no splicing
                end = start + len(parts[9])
                if strand=='0': end += extension
                else: start -= extension
            else: #splicing
                N=re.findall(r'[0-9]+', cigar) #numbers
                end = start + int(N[0])
                if strand=='0': end += extension
                else: start -= extension
            bases+=(end-start)
            try:
                starts[chro][start]+=1
                ends[chro][end]+=1
            except:
                starts[chro]=defaultdict(int)
                ends[chro]=defaultdict(int)
                starts[chro][start]=1
                ends[chro][end]=1
            if len(N)==3: #spliced alignment
                start2 = start + int(N[0])+ int(N[1])
                end2 = start2 + int(N[2])
                starts[chro][start2]+=1
                ends[chro][end2]+=1
                
        f.close()
    return starts, ends, unaligned, aligned, forward, reverse, gsize, bases

def printwig(outputfile, starts, ends, chrofile, outratio):
    allkeys={}
    f=open(outputfile, 'w')
    #print >>f, "track type=bedGraph name="+os.path.basename(outputfile)
    oof=open(chrofile, 'w')
    for chro in ends.keys():
        oof.write(chro+'\t'+str(max(ends[chro])+100000)+'\n')
    oof.close()
    print starts.keys()
    for chro in starts.keys():
        allkeys[chro]=list(set(starts[chro].keys())|set(ends[chro].keys()))
        allkeys[chro].sort()
    
    totalbase=0
    chrolen=1
    genolen=0
    print 'all', allkeys.keys()
    for chro in allkeys.keys():
        print 'chro     ', chro
        # ch=chro.upper()
        # if ch.startswith('CHRO'):
            # ch=ch.replace('CHRO', '')
        # elif ch.startswith('CHR'):
            # ch=ch.replace('CHR', '')
        # elif ch.startswith('CH'):
            # ch=ch.replace('CH', '')
        # if chro not in ['1','2','3','4','5','6','7','8','9','10','11', '12',\
		      # '13', '14', '15','16','17','18','19','20','21','22','X','Y','MT']:
		    # continue
        # try: size=chrosize[ch]
        # except: 
            # size=10000000000000
            # #print 'warning: did ot find chrosize for', ch, 'using', size
        counter=0
        prevkey=1
        allkey=allkeys[chro]
        start=starts[chro]
        end=ends[chro]
        #if chro.startswith('chr'): 
        print >>f, 'variableStep  chrom='+chro
        #else: print >>f, 'variableStep  chrom=chr'+chro
        for i in range(len(allkey)):
            key=allkey[i]
            #if key >=size: break
            #print chro, prevkey, key, counter
            #if counter!=0: print >> f, prevkey, counter/outratio
            print >> f, prevkey, counter/outratio
            totalbase+=(key-prevkey)*counter
            prevkey=key
            chrolen=key
            if start.has_key(key): counter+=start[key]
            if end.has_key(key): counter-=end[key]
        genolen+=chrolen
        print 'chro=', chro, 'chrolen=', chrolen, 'genolen=', genolen
    f.close()
    depth=float(totalbase)/genolen
    return depth

if __name__ == '__main__':
    (strand, extension, normalize, args)= options()
    #try: chrosize, chrofile = readchrosize(organism)
    #except: chrofile='/mnt/cluster/xdeng/script/ChromSizes/'+
    
    extension=int(extension)

    outputfile = args[0]+'.wig'
    dir=os.path.dirname(os.path.abspath(outputfile))
    chrofile=dir+'/'+str(uuid.uuid4())+'.size'
    inputfiles = args
    starts, ends, unaligned, aligned, forward, reverse, gsize, bases = loadsam(inputfiles, strand, extension)
    totalreads=unaligned+aligned
    if int(normalize)==0: outratio=1
    else: outratio=float(totalreads)/int(normalize)
    print 'total reads', totalreads, 'outratio', outratio
    depth = printwig(outputfile, starts, ends, chrofile, outratio)

    # print 'sample', inputfiles[0], 'depth', depth, 'unaligned', unaligned, 'aligned reads', aligned,'aligned percentage', round(float(aligned)/(aligned+unaligned),3)*100,'%'
    cmd = '/mnt/cluster/bsidna1_local/ChIPseeqer-2.1/dist/SCRIPTS/wigToBigWig -clip '+outputfile+' '+chrofile+' '+outputfile+'.bw'
    print cmd
    os.system(cmd)
    #print 'forward reads', forward
    #print 'reverse reads', reverse
    # try: print 'aligned percentage', round(float(aligned)/(aligned+unaligned),3)*100,'%'
    # except: pass
