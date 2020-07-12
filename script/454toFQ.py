#!/usr/bin/env python
import sys
f=open(sys.argv[1], 'r')
f2=open(sys.argv[2], 'r')
of=open(sys.argv[3],'w')

def phred2ASCII(line2):
	rval=[]
	quals = line2.strip().split()
	for q in quals:
		rval.append(chr(int(q) +33))
	return ''.join(rval)

seq=[]
qual=[]
id=''
i=0
for line1 in f:
	line2=f2.readline()
	if line1.startswith('>'):
		i+=1
		if id!='': 
			of.write(id+'\t 1\n'+''.join(seq)+'\n+\n'+''.join(qual)+'\n')
		id ='@'+str(i)
		seq=[]
		qual=[]
	else:
		seq.append(line1.strip())
		qual.append(phred2ASCII(line2))
#print id
of.write(id+'\t 1\n'+''.join(seq)+'\n+\n'+''.join(qual)+'\n')

f.close()
f2.close()
of.close()

