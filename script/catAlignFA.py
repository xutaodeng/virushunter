#!/usr/bin/python

import sys
import linecache
from operator import itemgetter, attrgetter

#print sys.argv[1]
e1=float(sys.argv[1])
e2=float(sys.argv[2])
f=open(sys.argv[3], 'r')
exaln=sys.argv[3]
exfa=sys.argv[4]
sortaln=open(sys.argv[5], 'w')
sortfa=open(sys.argv[6], 'w')

index=[]
i=0
offset=0
header=[]
for line in f:
	if line.strip().startswith('<td>') or line.strip().startswith('<tr>'):
		i+=1
		if line.strip().startswith('<tr>'):
			parts=line.strip()[4:].split('<td>')
		else:
			parts=line.strip().split('<td>')
		eval=float(parts[4])
		if e1<=eval and eval<=e2:
			index.append((i, eval))
	elif not line.strip().startswith('</table></html>'):
		header.append(line.strip())
		offset+=1
		print >>sortaln, line
	

index=sorted(index, key=itemgetter(1))
for (i, e) in index:
	sortaln.write(linecache.getline(exaln, i+offset).strip()+'\n')
	sortfa.write(linecache.getline(exfa, 2*i-1).strip()+'\n')
	sortfa.write(linecache.getline(exfa, 2*i).strip()+'\n')
sortaln.write('</table></html>\n')
f.close()
sortaln.close()
sortfa.close()
