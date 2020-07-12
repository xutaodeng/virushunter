#convert chipseqer output to SPP for IDR

import sys
file1, file2=sys.argv[1], sys.argv[2]
f=open(file1, 'r')
of=open(file2, 'w')

#########chipseeker format
# 1.Chromosome		: chromosome name
# 2.Start_Position		: the first (genomic) coordinate of the peak
# 3. End_Position		: the second (genomic) coordinate of the peak
# 4. Avg_p-value		: average log p-value of the nucleotides in the normalized peak region
# 5. Score			: score estimated as the average ChIP reads/length of peak (minus the average INPUT reads/length of peak - if INPUT is available)
# 6. Posmaxpeakheight   	: the position of the maximum height of the peak
# 7. Maxpeakheight		: the maximum peak height (in reads)
# 8. RelPosMaxPeakHeight(%)	: the relative position of the maximum height of the peak, e.g., 50% means the highest point is at the middle of the peak
# 9. Peak_Size		: the size of the peak (in bp)
# 10. Mid_point		: the middle position of the peak
# 11. Summit_dist_from_mid 	: the distance of the maximum height from the middle of the peak

############IDR format
# 1.chrom	 string	 Name of the chromosome
# 2.chromStart	 int	 The starting position of the feature in the chromosome. The first base in a chromosome is numbered 0.
# 3.chromEnd	 int	 The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the   feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
# 4.name	 string	 Name given to a region (preferably unique). Use '.' if no name is assigned
# 5.score	 int	 Indicates how dark the peak will be displayed in the browser (1-1000). If '0', the DCC will assign this based on signal value.         Ideally average signalValue per base spread between 100-1000.
# 6.strand	 char	 +/- to denote strand or orientation (whenever applicable). Use '.' if no orientation is assigned.
# 7.signalValue	 float	 Measurement of overall (usually, average) enrichment for the region.
# 8.pValue	 float	 Measurement of statistical signficance (-log10). Use -1 if no pValue is assigned.
# 9.qValue	 float	 Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned.
# 10.peak	 int	 Point-source called for this peak; 0-based offset from chromStart. Use -1 if no point-source called.

for line in f:
	parts=line.strip().split()
	chro, start, end, p, score, h1,h2,h3,h4,h5,h6=parts
	out=[chro,start,end,'.',score, '.',score,str(-float(p)),'.','-1']
	of.write('\t'.join(out)+'\n')
f.close()
of.close()
