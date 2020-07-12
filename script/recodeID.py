#!/usr/bin/env python
import sys
from collections import defaultdict
from operator import itemgetter, attrgetter

f=open(sys.argv[1], 'r')#fq dup file
of=open(sys.argv[2],'w') #sequence.txt with recoded read ID
label=sys.argv[3] #library label
try: pair_end =sys.argv[4] #pair 1 or 2
except: pair_end ='1'

i=0 #this is to be consistent with fq2faID.py
for line in f:
	i+=1
	if i%4==1:
		lineno=i/4
		#seqid='@s'+str(lineno)+' '+pair_end
		seqid='@s'+str(lineno)+'_'+pair_end+'_'+label
		print >>of, seqid
	else:
		print >>of, line.strip()
f.close()
of.close()