#!/usr/bin/env python
import sys
#samtools view -h -o out.sam  UNCID_358873.6e3fb6a1-d526-43b3-85f3-d995b5173fe5.110215_UNC3-RDR300156_00072_FC_62J8HAAXX.2.trimmed.annotated.translated_to_genomic.bam
f=open(sys.argv[1], 'r')
of=open(sys.argv[2], 'w')
unmap=False
try: 
	if sys.argv[3]=='unmap': unmap=True
except: pass
for line in f:
	if line.startswith('@'):continue
	parts=line.strip().split()
	id, ref, seq, qual=parts[0],parts[2],parts[9],parts[10]
	if ref!='*' and unmap: continue
	else: of.write('\n'.join(['@'+id, seq, '+', qual])+'\n')
f.close()
of.close()
