#!/usr/bin/env python
import sys
csvs=sys.argv[1:-2]

of1=open(sys.argv[-1],'w') #normalize
of2=open(sys.argv[-2], 'w') #read count
# header='miRNA'+'\t'+'\t'.join(csvs)
# of1.write(header+'\n')
# of2.write(header+'\n')

i=0

for csv in csvs:
	i+=1
	f=open(csv, 'r')
	f.readline() #header
	miRNAs,out1, out2=[],[],[]
	for line in f:
		parts=line.strip().split()
		miRNAs.append(parts[0])  #first file output miRNA 
		out1.append(parts[-1])
		out2.append(parts[-2])
	f.close()
	if i==1: 
		of1.write('Sample,'+','.join(miRNAs)+'\n')
		of2.write('Sample,'+','.join(miRNAs)+'\n')
	of1.write(csv[29:-4]+','+','.join(out1)+'\n')
	of2.write(csv[29:-4]+','+','.join(out2)+'\n')
of1.close()
of2.close()
