#!/usr/bin/env python
#  wig_peak_step1.py
# source wig2w2.sh
# source run.sh
#

import os
import os.path
import sys
from collections import defaultdict
from operator import itemgetter
import os
cwd = os.getcwd()

of1=open(cwd+'/wig2w2.sh', 'w')
of2=open(cwd+'/run.sh', 'w')

servers=['bsidna2','bsidna3', 'bsidna4','bsidna5','bsidna6','bsidna7','bsidna8','bsidna9','bsidna10','bsidna11',\
'bsidna12','bsidna13','bsidna14','bsidna16','bsidna17','bsidna18','bsidna19','bsidna20','bsidna22', \
'bsidna23', 'bsidna25',  'bsidna27','bsidna28','bsidna29', 'bsidna26']

dir='/mnt/cluster2/deng2/Kat/katerina/'

#2 peak files:
# allpeaks=['BCL6_5527_intron.inter.dist_hg38.txt',
# 'BCL6_5527_promoter_hg38.txt']

allpeaks=[
#'LSD1_2318_BCL6.nt-r1-r2-_1222.txt'
'LSD1_4300_promoter_hg38_1490.txt',
'LSD1_4300_intro.inter.dis_hg38_2346.txt'

# 'BCL6_5527_intron.inter.dist_hg38.txt',
# 'BCL6_5527_promoter_hg38.txt',
# 'LSD1_4300_intro.inter.dis_hg38_2346.txt',
# 'LSD1_4300_promoter_hg38_1490.txt',
# 'LSD1_2318_intro.inter.dis_hg38.txt',
# 'LSD1_2318_promoter_hg38.txt'
]


#8 wig files:
allwigs=[
'SUD4_H3K4me1_si_fc.wig',
'SUD4_H3K4me1_NT_fc.wig',
'HBL1_H3K4me1_si_fc.wig',
'HBL1_H3K4me1_NT_fc.wig',

'SUD4_BCL6_si_fc.wig',
'SUD4_BCL6_NT_fc.wig',
'SUD4_LSD1_si_fc.wig',
'SUD4_LSD1_NT_fc.wig',
'HBL1_BCL6_si_fc.wig',
'HBL1_BCL6_NT_fc.wig',
'HBL1_LSD1_si_fc.wig',
'HBL1_LSD1_NT_fc.wig'


# 'Sample_LSD1_NT_fc.wig2',
# 'Sample_LSD1_si_fc.wig2',
# 'Sample_H3K4me1_NT_fc.wig2',
# 'Sample_H3K4me1_si_fc.wig2',
# 'Sample_BCL6_NT_fc.wig2',
# 'Sample_BCL6_si_fc.wig2',
# 'Sample_H3K27ac_siCont_fc.wig2',
# 'Sample_H3K27ac_siBCL6_fc.wig2',
# 'Sample_H3K4me2_NT_fc.wig2',
# 'Sample_H3K4me2_si_fc.wig2'
]


allpeaks=[dir+pk for pk in allpeaks]

allwigs=[dir+pk for pk in allwigs]


cmd='python /mnt/cluster/xdeng/script/wig_peak_step2.py '

cmd2='python /mnt/cluster/xdeng/script/bedGraph2wig.py '#Sample_BCL6_NT_fc.wig Sample_BCL6_NT_fc.wig2 

i=0
for a in allwigs:
	i+=1
	of1.write('ssh '+servers[i%len(servers)]+' '+cmd2+' '+a+' '+a+'2 &\n')
of1.close()
i=0
for a in allwigs:
	for b in allpeaks:
		i+=1
		of2.write('ssh '+servers[i%len(servers)]+' '+cmd+' '+a+'2 '+b+' '+a+'_'+os.path.basename(b)+'.plot &\n')
of2.close()
#sys.exit(1)
