from collections import defaultdict
nuc =set(['A','C','G','T'])
def ti_tv(ref, alt):
	if not ((ref in nuc) and (alt in nuc)): return None
	if (ref=='A' and alt=='G') or (ref=='G' and alt=='A') or \
	(ref=='C' and alt=='T') or (ref=='T' and alt=='C'): return 'ti'
	else: return 'tv'
counts=defaultdict(int)
counts2=defaultdict(int)
f = open('snpEff.filtered.vcf', 'r')
i=0
for line in f:
	i+=1
	if i%100000 ==0: print i
	if line.startswith ('##'): continue
	if line.startswith ('#'): header = line; continue
	parts = line.strip().split()
	CHROM, POS,ID,REF,ALT,QUAL,FILTER,INFO=parts[0:8]
	rval  = ti_tv(REF, ALT)
	if rval == None: continue
	counts[rval]+=1
	if 'SYNONYMOUS' in INFO.upper(): counts2[rval]+=1

print 'genome: TI', counts['ti'], 'TV', counts['tv'], 'ti_tv_ratio', float(counts['ti'])/counts['tv']
print 'coding region: TI', counts2['ti'], 'TV', counts2['tv'], 'ti_tv_ratio', float(counts2['ti'])/counts2['tv']