import os
from os.path import join, getsize
f=open('SeqCap_EZ_Exome_v3_capture.bed','r')
of=open('SeqCap_EZ_Exome_v3_capture_fix.bed','w')
f.readline()
for line in f:
	if line.strip().startswith('chr'):
		of.write(line.strip()[3:]+'\n')
f.close()
of.close()
