#!/usr/bin/env python
import sys
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def PrintXML(subfile, of, a,b):
	nline=file_len(subfile)
	print subfile, nline
	f=open(subfile, 'r')
	i=0
	for line in f:
		i+=1
		if i>=a and i<=nline-b:
			print >>of, line.rstrip()
	f.close()

def mergeXML(filename,n):
	of=open(filename, 'w')
	for i in xrange(n):
		parts=filename.rsplit('.')
		subfile=parts[0]+'_'+str(i)+'.'+parts[1]
		if i==0:
			PrintXML(subfile, of, 1,2)
		elif i!=n-1:
			PrintXML(subfile, of, 21,2)
		else:
			PrintXML(subfile, of, 21,0)
	of.close()



infile=sys.argv[1]
n=int(sys.argv[2])
mergeXML(infile, n)

