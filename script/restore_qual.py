#!/usr/bin/env python
#clean header of fa files, only keep the first element
import sys
f=open(sys.argv[1], 'r') #fa file with indexed ID #blast.sig.fa files
f2=open(sys.argv[2], 'r') #fa file with complete ID
of=open(sys.argv[3], 'w') #fq file as output
of2=open(sys.argv[4], 'w') #fq file as output
i=0
allindex=set()
for line in f:
	i+=1
	if line.strip().startswith('>'):
		index=int(line.strip()[1:].strip()[1:])
		allindex.add(index)
i=0
for line in f2:
	i+=1
	if line.strip().startswith('>'):
		started=False
		if i in allindex:
			id, qid, qual=line.strip()[1:].strip().split()
			started=True
	elif started:
		seq=line.strip()
		print >>of, '\n'.join([id, seq, qid, qual])
		print >>of2, '\n'.join(['>'+id, seq])
f.close()
f2.close()
of.close()
of2.close()
#restore_qual.py blastn.sig.fa s_collapse.fa blastn.sig.fq