#!/usr/bin/env python
from collections import defaultdict
import os.path
import os
import string
import sys
import itertools

thread='30'
wd = os.path.abspath(os.path.dirname('.')).replace('san2', 'cluster2')

STAR='/mnt/cluster/xdeng/tools/STAR-STAR_2.4.0k/bin/Linux_x86_64/STAR'
samtools='/mnt/cluster/xdeng/tools/samtools-1.3/samtools'
genomeDir=wd+'/hg19_1/'
ref='/mnt/cluster/xdeng/gatk2.8Resource/human_g1k_v37.fasta'

# servers=[ 'bsidna2', 'bsidna3', 'bsidna4', \
 # 'bsidna5', 'bsidna6','bsidna7','bsidna8', \
  # 'bsidna9', 'bsidna10','bsidna11', 'bsidna12','bsidna13', \
  # 'bsidna14', 'bsidna15', 'bsidna16', 'bsidna17', \
# 'bsidna18', 'bsidna19', 'bsidna20', \
 # 'bsidna22', 'bsidna23', 'bsidna24', 'bsidna25', 'bsidna26', \
# 'bsidna27','bsidna28','bsidna29'
# ]
#servers=['bsidna32', 'bsidna33', 'bsidna34', 'bsidna35', 'bsidna36']
servers=['bsidna32', 'bsidna33',  'bsidna35']
nservers=len(servers)

def readSeeds_Single(seedfile='fastq/samples.txt'):
	seeds=defaultdict(list)
	f=open(seedfile, 'r')
	for line in f:
		key=line.strip()#split('_')[0]
		ss=line.strip()
		if line.strip().endswith('fasta'): ss = ss.replace('fasta', 'fastq')
		elif line.strip().endswith('sam'): ss = ss.replace('sam', 'fastq')
		elif line.strip().endswith('bam'): ss = ss.replace('bam', 'fastq')
		seeds[key].append(wd+'/fastq/'+ss)
	f.close()
	return seeds
########################################
# STAR
# #########################################
def STARrun1():
	#only run the first time to generate genome index
	cmd =STAR+' --runMode genomeGenerate --genomeDir '+genomeDir+' --genomeFastaFiles '+ref+' --runThreadN '+thread
	return cmd

def STARrun2(serverTag):
	cmd2= 'mkdir '+runDir1
	cmd3= 'cd '+runDir1
	cmd4= STAR+' --genomeDir '+genomeDir+' --readFilesIn '+read1+' '+read2+' --readFilesCommand zcat --runThreadN '+thread
	cmd5= 'mkdir '+genomeDir2
	cmd6=STAR+' --runMode genomeGenerate --genomeDir '+genomeDir2+' --genomeFastaFiles '+ref+' --sjdbFileChrStartEnd '+runDir1+'SJ.out.tab --sjdbOverhang 75 --runThreadN '+thread
	cmd7=' mkdir '+runDir2
	cmd8=' cd '+runDir2
	cmd9=STAR+' --genomeDir '+genomeDir2+' --readFilesIn '+read1+' '+read2+' --readFilesCommand zcat --runThreadN '+thread
	cmd=serverTag+ ' "'+' && '.join([cmd2, cmd3, cmd4, cmd5, cmd6, cmd7, cmd8, cmd9])+'" &\n'
	return cmd

def extract(serverTag,sample, ch, start, end):
	s4=runDir2+'s4.bam'
	s5=runDir2+sample+'extract1.sam'
	cmd = serverTag+' /mnt/cluster/xdeng/tools/samtools-1.3/samtools view '+s4+ ' "'+ch+':'+start+'-'+end+'" > '+s5+' &\n'
	return cmd

########################################
# GATK
#########################################
def GATKrun(serverTag):
	s1=runDir2+'Aligned.out.sam'
	s2=runDir2+'s2.bam'
	s3=runDir2+'s3.bam'
	s4=runDir2+'s4.bam'
	s5=runDir2+'s5.bam'
	s6=runDir2+'s6.bam'
	s7=runDir2+'s7.bam'
	out1=runDir2+'output1.vcf'
	out2=runDir2+'output2.vcf'
	out3=runDir2+'output3.vcf'
	outhtml=runDir2+'output3.html'
	recal=runDir2+'recal.table'
	dupfile=runDir2+'output.metrics'

	picard='/mnt/cluster/xdeng/tools/picard-tools-1.119/'
	gatk='/mnt/cluster/xdeng/tools/GenomeAnalysisTK-3.3-0/'
	dbSNP='/mnt/cluster/xdeng/gatk2.8Resource/dbsnp_138.b37.vcf'
	knownindel='/mnt/cluster/xdeng/gatk2.8Resource/Mills_and_1000G_gold_standard.indels.b37.vcf'

	option='-Djava.io.tmpdir=/mnt/cluster2/tmp  -Xms1024m -Xmx4096m '
	cmd1='java '+option +'-jar '+picard+'SortSam.jar I='+s1+' O='+s2+' SO=coordinate VALIDATION_STRINGENCY=SILENT'
	cmd2='java  '+option +'-jar '+picard+'AddOrReplaceReadGroups.jar I='+s2+' O='+s3+' SO=coordinate RGID=id RGLB=solexa-123 RGPL=ILLUMINA RGPU=AXL2342 RGSM='+sample
	cmd3='java '+option +'-jar '+picard+'MarkDuplicates.jar I='+s3+' O='+s4+'  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M='+dupfile
	cmd4=samtools +' index '+s4
	cmd5='java '+option +'-jar '+gatk+'GenomeAnalysisTK.jar -T SplitNCigarReads -R '+ref+' -I '+s4+' -o '+s5+' -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS'
	cmd6='java '+option +'-jar '+gatk+'GenomeAnalysisTK.jar -T BaseRecalibrator -nct '+thread+' -R '+ref+' -I '+s5+' -knownSites '+dbSNP+' -knownSites '+knownindel+' -o '+recal
	cmd7='java '+option +'-jar '+gatk+'GenomeAnalysisTK.jar -T PrintReads -nct '+thread+' -R '+ref+' -I '+s5+' -BQSR '+recal+' -o '+s6
	cmd8='java '+option +'-jar '+gatk+'GenomeAnalysisTK.jar -T HaplotypeCaller -R '+ref+' -I '+s6+' --dbsnp '+dbSNP+' -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o '+out1
	#cmd9='java '+option +'-jar '+gatk+'GenomeAnalysisTK.jar -T VariantFiltration -R '+ref+' -V '+out1+' -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o '+out2
	cmd10='java '+option +'-jar /mnt/cluster/tools/snpEff/snpEff.jar eff -c /mnt/cluster/tools/snpEff/snpEff.config -s '+outhtml+' -v GRCh37.74 '+out1+' > '+out3
	cmd=serverTag+ ' \''+' && '.join([cmd1, cmd2, cmd3, cmd4, cmd5, cmd6, cmd7, cmd8, cmd10])+'\' &\n'
	#cmd=serverTag+ ' \''+' && '.join([cmd9, cmd10])+'\' &\n'
	return cmd

if __name__ == "__main__":
	pair=False
	seeds=readSeeds_Single()
	try: os.mkdir(wd+'/STAR/')
	except: pass
	try: os.mkdir(genomeDir)
	except: pass
	job=0
	f11=open('star_run1.sh', 'w')
	f12=open('star_run2.sh', 'w')
	f2=open('gatk_run.sh', 'w')
	f3=open('bam_extract.sh', 'w')
	i=0
	vcfs=[]
	for sample in seeds.keys():
		if i==0:
			cmd=STARrun1()
			f11.write(cmd)
		i+=1
		serverTag = 'ssh '+servers[job%nservers]
		read1=seeds[sample][0]
		if pair: 
			read2=seeds[sample][1]
		else: read2=''
		runDir1=wd+'/STAR/'+sample+'/'
		runDir2=wd+'/STAR/'+sample+'/2pass/'
		vcfs.append(runDir2+'output3.vcf')
		genomeDir2=wd+'/STAR/'+sample+'/hg19_2/'
		cmd1=STARrun2(serverTag)
		cmd2=GATKrun(serverTag)
		cmd3=extract(serverTag, sample, '3', '38182601', '38182681')
		job+=1
		f12.write(cmd1)
		f2.write(cmd2)
		f3.write(cmd3)
		if job%nservers==0: 
			f12.write('wait\n')
			f2.write('wait\n')
			#f3.write('wait\n')
	f2.write('VCFsum.py '+' '.join(vcfs)+'\n')
	#f2.write('VCFdiffSum.py '+' '.join(vcfs)+'\n')
	f11.close()
	f12.close()
	f2.close()
	f3.close()
