import sys
import random
script,i,o1,o2=sys.argv
infile=open(i,'r')
out1=open(o1,'w')
out2=open(o2,'w')

for i in infile:
	file=random.randint(1,2)
	if i.startswith('#'):
		out1.write(i)
		out2.write(i)
		continue
	if file==1:
		out1.write(i)
	else:
		out2.write(i)
infile.close()
out1.close()
out2.close()
