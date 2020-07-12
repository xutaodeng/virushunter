#!/usr/bin/env python

import os
from collections import defaultdict
import os.path
import sys
import math

def average(s): return sum(s) * 1.0 / len(s)
def cut(s, c): return sum([1 for x in s if x>=c])
def sd(s):
	avg = average(s)
	variance = map(lambda x: (x - avg)**2, s)
	standard_deviation = math.sqrt(average(variance))
	return standard_deviation

f=open(sys.argv[1],'r')
of=open(sys.argv[2],'w')
sample=sys.argv[2].split('/')[-2]
header = f.readline()
i, Ns, means = 0,0,[]
for line in f:
	parts = line.strip().split()
	ref, N, mean =parts[0:3]
	N, mean = float(N), float(mean)
	Ns+=N
	means.append(mean)
	i+=1
print >>of, 'sample\tcov_avg\tcov_sd\tntarget\tPer10X\tPer20X\tPer50X\tPer100X'
print >>of, sample+'\t'+str(int(average(means)))+'\t'+str(int(sd(means)))+'\t'+str(i)+ \
		'\t'+str(cut(means,10)*100/i)+'\t'+str(cut(means, 20)*100/i)+'\t'+str(cut(means, 50)*100/i)+'\t'+str(cut(means, 100)*100/i)
of.close()