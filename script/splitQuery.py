#!/usr/bin/env python
import sys, os, os.path
def countFasta(filename):
    count=0
    f = open(filename, 'r')
    for line in f:
        if line.strip().startswith('>'):
            count+=1
    f.close()
    return count 

def splitFasta(filename, nfiles, nseqs):
    # chunksize = int(float(nseqs)/nfiles+1)
    fs= []
    for i in xrange(nfiles):
        fname=filename+'_'+str(i)
        f = open(fname, 'w')
        fs.append(f)
    f = open(filename, 'r')
    handle =0
    j=0 #index in handle
    for line in f:
        if line.strip().startswith('>'):
            j+=1
        fs[j%nfiles].write(line)

    for i in xrange(nfiles):
        fs[i].close()
    f.close()
    

infile=sys.argv[1]
n=int(sys.argv[2])
nseqs=countFasta(infile)
print infile, 'num_seqs_to_blast = ', nseqs
splitFasta(infile, n, nseqs)

