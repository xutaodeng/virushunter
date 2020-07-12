#!/usr/bin/env python
import sys
import os
infile = sys.argv[1]
print infile

cmd = 'samtools view -bS -o '+infile+'.bam '+infile
print 'executing', cmd
os.system(cmd)

cmd = 'samtools sort '+infile+'.bam '+infile+'.sorted'
print 'executing', cmd
os.system(cmd)

cmd = 'samtools pileup -cv -f /home/xdeng/bowtie-0.11.3/indexes/Homo_sapiens.GRCh37.56.dna.toplevel.fa '+infile+'.sorted.bam > '+infile+'.pileup'
print 'executing', cmd
os.system(cmd)

cmd = 'samtools.pl varFilter '+infile+'.pileup | awk '+''' '$6>=20 && $8>=20' > ''' +infile+'.final.pileup'
print 'executing', cmd
os.system(cmd)
