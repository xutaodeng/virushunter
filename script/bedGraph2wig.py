#!/usr/bin/env python
from collections import defaultdict
import os.path
import os
import string
import sys
import itertools

# bedGraph2wig.py Sample_BCL6_NT_fc.wig Sample_BCL6_NT_fc.wig2 &
# bedGraph2wig.py Sample_BCL6_si_fc.wig Sample_BCL6_si_fc.wig2 &
# bedGraph2wig.py Sample_H3K27ac_siBCL6_fc.wig Sample_H3K27ac_siBCL6_fc.wig2 &
# bedGraph2wig.py Sample_H3K27ac_siCont_fc.wig Sample_H3K27ac_siCont_fc.wig2 &
# bedGraph2wig.py Sample_H3K4me1_NT_fc.wig Sample_H3K4me1_NT_fc.wig2 &
# bedGraph2wig.py Sample_H3K4me1_si_fc.wig Sample_H3K4me1_si_fc.wig2 &
# bedGraph2wig.py Sample_H3K4me2_NT_fc.wig Sample_H3K4me2_NT_fc.wig2 &
# bedGraph2wig.py Sample_H3K4me2_si_fc.wig Sample_H3K4me2_si_fc.wig2 &
# bedGraph2wig.py Sample_LSD1_NT_fc.wig Sample_LSD1_NT_fc.wig2 &
# bedGraph2wig.py Sample_LSD1_si_fc.wig Sample_LSD1_si_fc.wig2 &


f=open(sys.argv[1], 'r') # begGraph format
of = open(sys.argv[2], 'w') #wig fixsize format
span=10
  # variableStep  chrom=chrN  [span=windowSize]
  # chromStartA  dataValueA
  # chromStartB  dataValueB
  # ... etc ...  ... etc ...

chros=set([])
for line in f:
	if line.startswith('#'): continue
	parts=line.strip().split()
	if parts[3]==0 : continue
	chro, start, end, val = parts
	if chro not in chros:
		of.write('variableStep  chrom='+chro+'  span='+str(span)+'\n')
		chros.add(chro)
	start, end, val = int(start), int(end), float(val)
	pos = start
	while pos <= end:
		if pos==0: pos=1
		of.write(str(pos)+'\t'+str(val)+'\n')
		pos+=span
f.close()
of.close()