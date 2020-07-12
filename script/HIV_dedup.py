#!/usr/bin/env python

from os import walk
from collections import defaultdict
from operator import itemgetter, attrgetter
import sys

#nn=200
nn=None #full length

def dedup(infile, outfile, rate_thres=0.01):
	ff = open(infile, 'r')
	of =open(outfile, 'w')
	reads2=defaultdict(list)
	seq=[]
	name=''
	i=0
	for line in ff:
		if line.strip().startswith('>'):
			i+=1
			#if i%1000==0:print i, len(reads2)
			if name!='':
				#reads[''.join(seq)].append(name)
				frag = ''.join(seq)[1:nn]
				#if reads2.has_key(frag): continue
				reads2[frag].append(name)
			name=line.strip()
			seq=[]
		else:
			seq.append(line.strip())
	#reads[''.join(seq)].append(name)
	frag = ''.join(seq)[1:nn]
	reads2[frag].append(name)
	
	reads = reads2.items()
	freq_thres = int(i*rate_thres)
	if freq_thres<3: freq_thres=3
	rreads = [(read, names[0], len(names)) for (read, names) in reads if len(names)>=freq_thres ]
	print infile, 'n_reads', i, 'n_uniq', len(reads2), 'n_uniq_freq>=',freq_thres, len(rreads)
	rreads.sort(key=itemgetter(2), reverse=True)
	for read, name, freq in rreads:
		print >>of, name,freq,'\n'+read
	ff.close()
	of.close()

if __name__ == '__main__':
	infile = sys.argv[1]
	outfile =sys.argv[2] 
	rate_thres = float(sys.argv[3]) #default 0.01
	dedup(infile, outfile, rate_thres)

# f = []
# for (dirpath, dirnames, filenames) in walk('FastaOnly'):
    # f.extend(filenames)
    # break

# #nn=200
# nn=None
# for fi in f:
    # #if not fi.startswith('A-9'): continue

        

