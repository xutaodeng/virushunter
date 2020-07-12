#!/usr/bin/env python
from collections import defaultdict
import operator
import sys
import os
import os.path
from optparse import OptionParser

def options():
    print '''usage: %prog [--blast blastprogram ]
             blast against two fa files, you may limit query sequence length with -length'''

    args = sys.argv[2:]
    blast = sys.argv[1]
    length = 1
    buf ='' # the rest of args to be cascade to blast+
    subject =None
    out=None
    i=0
    while i < len(args):
        if args[i]== '-length': length = int(args[i+1]); i+=1
        elif args[i]=='-query': query = args[i+1];  i+=1
        elif args[i]=='-subject': subject = args[i+1]; i+=1
        elif args[i]=='-out': buf+=args[i]+' "'+args[i+1]+'" '; i+=1
        else: buf+=args[i]+' '
        i+=1
    return (blast, length, buf, query, subject)

def filterLength(infile, l):
	print 'length', l
	outfile = infile.rsplit('.', 1)[0]+str(l)+'_filter.fa'
	print 'converting', infile, 'to', outfile
	f = open(infile, 'r')
	of = open(outfile, 'w')
	i = 0
	seq=[]
	id=''
	a, b = 0, 0
	for line in f:
		if line.strip().startswith('>'):
			a+=1
			if id!='' and len (''.join(seq)) >= l:
				b+=1
				of.write('>'+id+'\n'+''.join(seq)+'\n')
			id = line.strip()[1:]
			seq=[]
		else:
			seq.append(line.strip())
	of.write('>'+id+'\n'+''.join(seq)+'\n')
	f.close()
	of.close()
	print a, b
	return outfile
	
if __name__ == '__main__':
	(blast, length, buf, query, subject)= options()
	#subject, query, output = args
	if subject!=None:
		subject = os.path.abspath(subject)
		buf += ' -subject "'+subject+ '" '
	if int(length)>1: query = filterLength(query, int(length))
	query = os.path.abspath(query)
	buf += ' -query "'+query+ '" '
	cmd = blast+' '+buf
	print 'running:', cmd
	os.system(cmd)