import gzip
from collections import defaultdict

count=defaultdict(int)
count2=defaultdict(set)
count3=defaultdict(set)

f=open('species_index.txt', 'r')
for line in f:
	ind, name = line.strip().split()
	cat, clas, fam, spec = name.split('_', 3)
	count[cat]+=1
	count2[cat].add(clas)
	count3[cat].add(fam)

f.close()
for cat in count.keys():
	print cat, 'spec', count[cat], 'fam', len(count3[cat]), 'class', len(count2[cat])
