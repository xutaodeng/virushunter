#!/usr/bin/env python
import os
from collections import defaultdict
import os.path
import sys

#rsync -ruv /source/directory/* username@domain.net:/destination/directory
#tail -n +2 SeqCap_EZ_Exome_v3_capture.bed > SeqCap_EZ_Exome_v3_capture_noheader.bed
#sed -i 's/chr//g' SeqCap_EZ_Exome_v3_capture_noheader.bed
#java -Xmx16g -jar /mnt/cluster/tools/BAMStats-1.25/BAMStats-1.25.jar -i s2.bam -o s2.summary -f /mnt/san2/xdeng/measles/ihg-client.ucsf.edu/seielstadm/SeqCap_EZ_Exome_v3_capture_fix.bed

#java  -Xms1024m -Xmx4096m -jar  snpEff.jar download GRCh37.64
#java  -Xms1024m -Xmx4096m -jar  snpEff.jar databases|grep 'GRCh37'
# http://www.broadinstitute.org/gatk/guide/topic?name=faqs 
# password: <blank>
# ftp gsapubftp-anonymous@ftp.broadinstitute.org
# mget dbsnp_138.b37.vcf.gz
# mget Mills_and_1000G_gold_standard.indels.b37.vqcf.gz
# mget 1000G_phase1.snps.high_confidence.b37.vcf.gz
# mget 1000G_omni2.5.b37.vcf.gz
# mget human_g1k_v37.fasta.gz
#java  -Xms1024m -Xmx4096m -jar  /mnt/cluster/tools/picard-tools-1.105/SamToFastq.jar INPUT=NA12892.bam FASTQ=1.fq SECOND_END_FASTQ=2.fq VALIDATION_STRINGENCY=SILENT

# bwa index -a bwtsw human_g1k_v37.fasta
#java  -Xms1024m -Xmx4096m -jar  /mnt/cluster/tools/picard-tools-1.105/CreateSequenceDictionary.jar R=human_g1k_v37.fasta O=human_g1k_v37.dict
#samtools faidx human_g1k_v37.fasta

ref = '/mnt/cluster2/gatk_resource2.8/human_g1k_v37.fasta'
picard = '/mnt/cluster/tools/picard-tools-1.105/'
vcfConcat='/mnt/cluster/tools/vcftools_0.1.13/bin/vcf-concat '
vcfMerge='/mnt/cluster/tools/vcftools_0.1.13/bin/vcf-merge '
gatk= '/mnt/san/cluster/tools/GenomeAnalysisTK-3.1-1/'
snpeff='/mnt/san/cluster/tools/snpEff/'
dbSNP='/mnt/cluster2/gatk_resource2.8/dbsnp_138.b37.vcf'
knownindel='/mnt/cluster2/gatk_resource2.8/Mills_and_1000G_gold_standard.indels.b37.vcf'
hapmap = '/mnt/cluster2/gatk_resource2.8/hapmap_3.3.b37.vcf'
omni = '/mnt/cluster2/gatk_resource2.8/1000G_omni2.5.b37.vcf'
onekg = '/mnt/cluster2/gatk_resource2.8/1000G_phase1.snps.high_confidence.b37.vcf'
bed ='/mnt/san2/xdeng/measles/SeqCap_EZ_Exome_v3_capture_fix.bed'
thread = '16'
servers=[ 'bsidna3','bsidna5','bsidna6','bsidna7','bsidna8','bsidna9','bsidna10',\
'bsidna11','bsidna17','bsidna18','bsidna19',\
'bsidna22','bsidna23','bsidna24','bsidna25','bsidna27','bsidna29', 'bsidna30']
option='-Djava.io.tmpdir=/mnt/cluster2/tmp  -Xms1024m -Xmx4096m '

bigservers=['bsidna33','bsidna34','bsidna35']

nservers=len(servers)
nbigservers=len(bigservers)
wd = os.path.abspath(os.path.dirname('.'))
print 'hello'


def bam2bigwig(sams):
	# f0 = open('samindex.sh', 'w')
	# f = open('bam2sam.sh', 'w')
	f2 = open('sam2wig.sh', 'w')
	#f3 = open('wig2bigwig.sh', 'w')
	job=0
	#if region is '', whole genome will be processed
	for sam in sams:
		serverTag = 'ssh '+servers[job%nservers]
		# base=wd+'/'+key+'_t_out/accepted_hits'
		# f0.write(serverTag+' "/mnt/cluster/xdeng/tools/samtools-1.3/samtools index '+base+'.bam" &\n') 
		# f.write(serverTag+' "/mnt/cluster/xdeng/tools/samtools-1.3/samtools view -h -o '+base+'_'+key+'.sam '+base+'.bam "'+region+'" "&\n')
		# #f.write(serverTag+' "/mnt/cluster/xdeng/tools/samtools-1.3/samtools view -h -o '+base+'.sam '+base+'.bam "&\n')
		f2.write(serverTag+' "/mnt/cluster/xdeng/script/sam2wig.py -n 1000000  '+sam+' "&\n') #no normalize
		# f3.write(serverTag+' "/mnt/san/cluster/bsidna1_local/ChIPseeqer-2.1/dist/SCRIPTS/wigToBigWig -clip '+ \
				# base+'.sam.wig '+chrofile+' '+base+'.sam.wig.bw"&\n')
		job+=1
		if job%nservers==0:  
			# f.write('wait\n')
			f2.write('wait\n')
			# f0.write('wait\n')
	# f.close()
	f2.close()
	# f0.close()
	#f3.close()
	
	
def samsort(sams):
	f1 = open('sam2bam.sh', 'w')
	f = open('bamsort.sh', 'w')
	f2 = open('split.sh', 'w')
	job=0
	i=0
	bams=[]
	for infile in sams:
		sortoutfile = infile+'.sort'
		outfile = infile+'.bam'
		#bams.append(infile+'.sort.bam')
		bams.append(infile+'.sort.bam')
		serverTag = 'ssh '+servers[job%nservers]
		f1.write(serverTag+' /mnt/cluster/xdeng/tools/samtools-1.2/samtools view -b -S -T '+ref+' '+infile+' > '+outfile+' &\n')
		f.write(serverTag+' /mnt/cluster/xdeng/tools/samtools-1.2/samtools sort  '+outfile+' '+sortoutfile+' &\n')
		f2.write(serverTag+' /mnt/cluster/tools/bamtools/build/mnt/cluster/tools/bamtools/build/bin/bamtools split -in '+sortoutfile+'.bam -reference &\n')
		job+=1
		if job%nservers==0: 
			f.write('wait\n')
			f1.write('wait\n')
		i+=1
	f.close()
	f1.close()
	f2.close()
	return sams


def mpileup(sams):
	f = open('pileup.sh', 'w')
	f2 = open('varscan.sh', 'w')
	f3 = open('snpeff.sh', 'w')
	f4 = open('vcf_merge.sh', 'w')
	f5 = open('vcf_sum.sh', 'w')
	f6 = open('dbsnp.sh', 'w')
	vcfs=[]
	indels=[]
	try: os.mkdir('pileup')
	except: pass
	job=0
	for chro in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10','11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT']:
	#for chro in ['4', '5', '6', '7', '8', '9', '10','11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT']:
		subs=[]
		for sam in sams:
			sub=sam+'.sort.REF_'+chro+'.bam'
			subs.append(sub)
		bigserverTag = 'ssh '+bigservers[job%nbigservers]
		serverTag = 'ssh '+servers[job%nservers]
		#f2 = open('consensus.sh', 'w')
		f.write(serverTag+' /mnt/cluster/xdeng/tools/samtools-1.2/samtools mpileup  -d 1000000 -C50 -Df '+ref+' ')
		f.write(' '.join(subs)+' > '+os.getcwd()+'/pileup/'+chro+'.pileup &\n')
		job+=1
		if job%nservers==0: 
			f.write('wait\n')
		snpvcf=os.getcwd()+'/pileup/'+chro+'.snp.vcf '
		indelvcf=os.getcwd()+'/pileup/'+chro+'.indel.vcf '
		snpvcf2=os.getcwd()+'/pileup/'+chro+'.snp2.vcf '
		indelvcf2=os.getcwd()+'/pileup/'+chro+'.indel2.vcf '
		snphtml=os.getcwd()+'/pileup/'+chro+'.snp.html '
		indelhtml=os.getcwd()+'/pileup/'+chro+'.indel.html '
		#f2 = open('consensus.sh', 'w')
		f2.write(serverTag+' java -jar /mnt/cluster/tools/VarScan.v2.3.9.jar mpileup2snp ')
		f2.write(os.getcwd()+'/pileup/'+chro+'.pileup --p-value 0.05 --min-coverage 5 --min-avg-qual 15 --min-var-freq 0.01 --output-vcf 1 ')
		f2.write(' > '+snpvcf+ ' &\n')
		f2.write(serverTag+' java -jar /mnt/cluster/tools/VarScan.v2.3.9.jar mpileup2indel ')
		f2.write(os.getcwd()+'/pileup/'+chro+'.pileup --p-value 0.05 --min-coverage 5 --min-avg-qual 15 --min-var-freq 0.01 --output-vcf 1 ')
		f2.write(' > '+indelvcf+' &\n')
		f3.write(serverTag+' java '+option +'-jar /mnt/cluster/tools/snpEff/snpEff.jar eff -c /mnt/cluster/tools/snpEff/snpEff.config -s '+snphtml+' -v GRCh37.74 '+snpvcf+' > '+snpvcf2+' &\n')
		f3.write(serverTag+' java '+option +'-jar /mnt/cluster/tools/snpEff/snpEff.jar eff -c /mnt/cluster/tools/snpEff/snpEff.config -s '+indelhtml+' -v GRCh37.74 '+indelvcf+' > '+indelvcf2+' &\n')
		vcfs.append(snpvcf2)
		indels.append(indelvcf2)
	f4.write('export PERL5LIB=/mnt/cluster/tools/vcftools_0.1.13/perl:$PERL5LIB\n')
	f4.write(vcfConcat+ ' '.join(vcfs)+' > all_snp.vcf &\n')
	f4.write(vcfConcat+ ' '.join(indels)+' > all_indel.vcf &\n')

	f6.write('VCF_SNPadd.py '+dbSNP+' all_snp.vcf  all_rs_snp.vcf &\n')
	f6.write('VCF_SNPadd.py '+dbSNP+' all_indel.vcf  all_rs_indel.vcf &\n')
	f5.write('VCFsumVarscan.py '+' '.join(samples.keys())+' \n')
	f.close()
	f2.close()
	f3.close()
	f4.close()
	f5.close()
	f6.close()
	bam2bigwig(sams)

# BSI\306307@bsidna1:/mnt/cluster2/deng2/JAMES_WES_human> VCFmerge2.py vcfsummarySNP.txt bach.2/vcfsummarySNP.txt vcfsummarySNP_combine.txt
# num mut 421684
# num samp 20
# num mut 1366617
# num samp 8
# total variants 1742326
# common variants 45975
# BSI\306307@bsidna1:/mnt/cluster2/deng2/JAMES_WES_human> VCFmerge2.py vcfsummaryIndel.txt bach.2/vcfsummaryIndel.txt vcfsummaryIndel_combine.txt
# num mut 19502
# num samp 20
# num mut 49934
# num samp 8
# total variants 63750
# common variants 5686

def readlanes():
	of=open('cat.sh','w')
	f=open('fastq/samples.txt', 'r')
	lanesR1=defaultdict(list)
	lanesR2=defaultdict(list)
	for line in f:
		if '_R1_' in line:
			key= line.strip().rsplit('_',3)[0]
			lanesR1[key].append(line.strip())
		if '_R2_' in line:
			key= line.strip().rsplit('_',3)[0]
			lanesR2[key].append(line.strip())
	i=0
	for key in lanesR1.keys():
		print key, lanesR1[key]
		print key, lanesR2[key]
		r1=lanesR1[key]
		r2=lanesR2[key]
		r1.sort()
		r2.sort()
		assert len(r1)==len(r2)
		server='ssh '+servers[i%nservers]+' '
		cmd1 = server+' "cd '+wd+'/fastq/'+' && zcat '+' '.join(r1)+' >'+key+'_1.fq" &'
		of.write(cmd1+'\n')
		i+=1
		# if i%nservers==0: 
			# of.write('wait\n')
		server='ssh '+servers[i%nservers]+' '
		cmd2 = server+' "cd '+wd+'/fastq/'+' && zcat '+' '.join(r2)+' > '+key+'_2.fq" &'
		of.write(cmd2+'\n')
		i+=1
		# if i%nservers==0: 
			# of.write('wait\n')
	f.close()
	of.close()

def readfq():
	samples=defaultdict(list)
	f=open('fastq/samples.txt', 'r')
	for line in f:
		key= line.strip().rsplit('_',4)[0]
		if not samples.has_key(key):
			samples[key]={1:[], 2:[]}
		if '_R1_' in line:
			samples[key][1].append(line.strip())
		if '_R2_' in line:
			samples[key][2].append(line.strip())

			# of.write('wait\n')
	f.close()
	for key in samples.keys():
		print key, samples[key]
	print 'numSamples', len(samples)
	return samples
#readlanes()
samples=readfq()

#lanes= os.listdir('.') 
# lanes = ["Project_seielstadm-MA34","Project_seielstadm-MA35","Project_seielstadm-MA36","Project_seielstadm-MA37",\
# "Project_seielstadm-MA38","Project_seielstadm-MA39","Project_seielstadm-MA40", "Project_seielstadm-MA41",
# "Project_seielstadm-MA42", \
# "Project_seielstadm-MA43", \
# "Project_seielstadm-MA44", \
# "Project_seielstadm-MA45", \
# "Project_seielstadm-MA46", \
# "Project_seielstadm-MA47", \
# "Project_seielstadm-MA48", \
# "Project_seielstadm-MA49", \
# "Project_seielstadm-MA50", \
# "Project_seielstadm-MA51", \
# "Project_seielstadm-MA52", \
# "Project_seielstadm-MA53", \
# "Project_seielstadm-MA54", \
# "Project_seielstadm-MA55" \
# ]

# lanesample = defaultdict(list)
# i=0
# for lane in lanes:
	# if os.path.isdir(lane):
		# samples= os.listdir(lane)
		# print lane
		# for sample in samples:
			# lanesample[lane].append(sample)
			# print '\t', sample
			# sname=sample.split('_')[1]
			# splits= os.listdir(lane+'/'+sample)
			# r1,r2=[],[]
			# for split in splits:
				# if 'R1' in split and split.endswith('.gz'):
					# print '\tRead1', split
					# r1.append(lane+'/'+sample+'/'+split)
				# if 'R2' in split and split.endswith('.gz'):
					# print '\tRead2', split
					# r2.append(lane+'/'+sample+'/'+split)
			# r1.sort()
			# r2.sort()
			# assert len(r1)==len(r2)
			# server='ssh '+servers[i%nservers]+' '
			# cmd1 = server+' "cd '+wd+' && zcat '+' '.join(r1)+' > '+lane+'/'+sample+'/1.fq" &'
			# f.write(cmd1+'\n')
			# i+=1
			# if i%nservers==0: 
				# f.write('wait\n')
			# server='ssh '+servers[i%nservers]+' '
			# cmd2 = server+' "cd '+wd+' && zcat '+' '.join(r2)+' > '+lane+'/'+sample+'/2.fq" &'
			# f.write(cmd2+'\n')
			# i+=1
			# if i%nservers==0: 
				# f.write('wait\n')
# f.close()

f=open('gatk_pipe_all.sh', 'w')
#dirs = ['NA12878', 'NA12891', 'NA12892']
#dirs= ['exampleBAM']
fam='input.fam' #first 6 columns of ped
i=0
dirs=[]
sams=[]
nbigservers=len(bigservers)
job=0
of=open('bowtie.sh', 'w')
for sample in samples.keys():
	fq1=[os.getcwd()+'/fastq/'+x for x in samples[sample][1]]
	fq2=[os.getcwd()+'/fastq/'+x for x in samples[sample][2]]
	bow='ssh '+bigservers[job%nbigservers]+' /mnt/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2 --quiet --local --no-hd --reorder -p 16 -x /mnt/cluster2/gatk_resource2.8/v37'+ \
	' -1 '+','.join(fq1)+' -2 '+','.join(fq2)+' -S '+os.getcwd()+'/'+sample+'.sam'
	sams.append(os.getcwd()+'/'+sample+'.sam')
	job+=1
	of.write(bow+' &\n')
	if job%nbigservers==0: 
		of.write('wait\n')
of.close()
bams=samsort(sams)
mpileup(bams)
ff=open('pipeline.sh', 'w')

cmd=['source bowtie.sh', 'source sam2bam.sh', 'source bamsort.sh','source split.sh','source pileup.sh','source varscan.sh',\
'source snpeff.sh','source vcf_merge.sh','source dbsnp.sh','source vcf_sum.sh']
ff.write('\nwait\n'.join(cmd)+'\n')
ff.close()





sys.exit()
	
for sample in samples.keys():
	dir =sample
	fq1=[os.getcwd()+'/fastq/'+x for x in samples[sample][1]]
	fq2=[os.getcwd()+'/fastq/'+x for x in samples[sample][2]]
	try: os.mkdir(dir)
	except: pass
	print dir
	dirs.append(dir)
	server='ssh '+servers[i%nservers]+' '
	ff=open(dir+'/gatk_pipe.sh', 'w')
	f.write('cd '+dir+'\n')
	wd=os.getcwd()+'/'+dir+'/'
	f.write(server+' source '+wd+'gatk_pipe.sh &\n')
	f.write( 'cd ..\n')
	#inputbam=dir+'.bam'
	sampleID=dir#inputbam.split('.')[0]
	#sam2fq = 'java  -Xms1024m -Xmx8g -jar  '+picard+'SamToFastq.jar INPUT='+wd+inputbam+' FASTQ='+wd+'1.fq SECOND_END_FASTQ='+wd+'2.fq VALIDATION_STRINGENCY=SILENT'
	# bwa1 = '/mnt/cluster/tools/bwa-0.5.9/bwa aln -t '+thread+' '+ref+ ' '+fq1+' > '+wd+'1.sai'
	# err = "rc=$? \nif [[ $rc != 0 ]] ; \nthen \nexit $rc \nfi\n"
	# bwa2 = '/mnt/cluster/tools/bwa-0.5.9/bwa aln -t '+thread+' '+ref+ ' '+fq2+' > '+wd+'2.sai'
	# bwape='/mnt/cluster/tools/bwa-0.5.9/bwa sampe '+ref+' '+wd+'1.sai '+wd+'2.sai '+fq1+' '+fq2+' > '+wd+'s1.sam'
	bow=' /mnt/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2 --quiet --local --no-hd --reorder -p 16 -x /mnt/cluster2/gatk_resource2.8/v37'+ \
			' -1 '+' '.join(fq1)+' -2 '+' '.join(fq2)+' -S '+wd+'s1.sam'
	err = "rc=$? \nif [[ $rc != 0 ]] ; \nthen \nexit $rc \nfi\n"
	sortbam ='java -Djava.io.tmpdir=/mnt/cluster2/tmp  -Xms1024m -Xmx4096m -jar  '+picard+'SortSam.jar I='+wd+'s1.sam O='+wd+'s2.bam SO=coordinate VALIDATION_STRINGENCY=SILENT'
	sams.append(wd+'s1.sam')
	bamstat= 'java -Djava.io.tmpdir=/mnt/cluster2/tmp  -Xmx8g -jar /mnt/cluster/tools/BAMStats-1.25/BAMStats-1.25.jar -i '+wd+'s2.bam -o '+wd+'s2.summary -f '+bed
	bamstatpy = '/mnt/cluster2/xdeng/measles/bamstat.py '+wd+'s2.summary '+wd+'s2.stat'
	dup ='java -Djava.io.tmpdir=/mnt/cluster2/tmp  -Xms1024m -Xmx4096m -jar  '+picard+'MarkDuplicates.jar I='+wd+'s2.bam O='+wd+'s3.bam METRICS_FILE='+wd+'metric.dup VALIDATION_STRINGENCY=SILENT'
	rg ='java -Djava.io.tmpdir=/mnt/cluster2/tmp  -Xms1024m -Xmx4096m -jar  '+picard+'AddOrReplaceReadGroups.jar I='+wd+'s3.bam O='+wd+'s4.bam RGID=id RGLB=solexa-123 RGPL=illumina RGPU=AXL2342 RGSM='+sampleID+' RGCN=kaiser RGDT=2014 VALIDATION_STRINGENCY=SILENT'
	index ='samtools index '+wd+'s4.bam'
	realign1 = 'java -Djava.io.tmpdir=/mnt/cluster2/tmp  -Xms1024m -Xmx4096m -jar  '+gatk+'GenomeAnalysisTK.jar -T RealignerTargetCreator -R '+ref+' -I '+wd+'s4.bam -known '+knownindel+' -o '+wd+'s.intervals'
	realign2 = 'java -Djava.io.tmpdir=/mnt/cluster2/tmp  -Xms1024m -Xmx4096m -jar  '+gatk+'GenomeAnalysisTK.jar -T IndelRealigner -R '+ref+' -I '+wd+'s4.bam -known '+knownindel+' -targetIntervals '+wd+'s.intervals -o '+wd+'s5.bam'
	bqsr1='java -Djava.io.tmpdir=/mnt/cluster2/tmp  -Xms1024m -Xmx4096m -jar  '+gatk+'GenomeAnalysisTK.jar -T BaseRecalibrator -nct '+thread+' -R '+ref+' -I '+wd+'s5.bam -knownSites '+dbSNP+' -knownSites '+knownindel+' -o '+wd+'recal.table'
	bqsr2='java -Djava.io.tmpdir=/mnt/cluster2/tmp  -Xms1024m -Xmx4096m -jar  '+gatk+'GenomeAnalysisTK.jar -T PrintReads -nct '+thread+' -R '+ref+' -I '+wd+'s5.bam -BQSR '+wd+'recal.table -o '+wd+'s6.bam'
	reduce = 'java -Djava.io.tmpdir=/mnt/cluster2/tmp  -Xms1024m -Xmx4096m -jar  '+gatk+'GenomeAnalysisTK.jar -T ReduceReads -R '+ref+' -I '+wd+'s6.bam -o '+wd+'s7.bam'
	#cmds = [sam2fq, bwa1, bwa2, bwape, sortbam, dup, rg, index, realign1, realign2, bqsr1, bqsr2, reduce]
	#cmds = [bwa1, err, bwa2, err, bwape, err, sortbam, bamstat, bamstatpy, dup, rg, index, realign1, realign2, bqsr1, bqsr2, reduce]
	cmds = [bow, err, sortbam, bamstat, bamstatpy, dup, rg, index, realign1, realign2, bqsr1, bqsr2, reduce] #bowtie
	#cmds = [bamstat, bamstatpy]
	#cmds = [bamstatpy]
	#msgs = ['echo '+cmd for cmd in cmds]
	#ff.write('\n'.join([x for y in zip(msgs,cmds) for x in y])+'\n')
	ff.write('\n'.join(cmds)+'\n')
	i+=1
	ff.close()
	if i%nservers==0: 
		f.write('wait\n')
bams=samsort(sams)
mpileup(bams)

f.write('wait\n')
wd = os.path.abspath(os.path.dirname('.'))
bams=' '.join(['-I '+wd+'/'+dir+'/s7.bam' for dir in dirs])
stats=' '.join([dir+'/s2.stat' for dir in dirs])
allstats = 'cat '+stats+' > all.stat'

chros=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','MT']

i=0
vcfs=[]
for chro in chros:
	call = 'ssh '+servers[i]+ ' java -Djava.io.tmpdir=/mnt/cluster2/tmp  -Xms1024m -Xmx8g -jar  '+gatk+'GenomeAnalysisTK.jar -nct '+thread+ \
		' -L '+chro+'  -T UnifiedGenotyper -glm BOTH -R '+ref+' '+bams+ \
		' -o '+wd+'/'+chro+'_raw.SNPs.vcf -stand_call_conf 30 -stand_emit_conf 10'
	vcfs.append(chro+'_raw.SNPs.vcf')
	#f2.write(call+' &\n')#
	i+=1
mergeVCF = 'java -Djava.io.tmpdir=/mnt/cluster2/tmp  -Xms1024m -Xmx4096m -jar  '+picard+'MergeVcfs.jar INPUT='+' INPUT='.join(vcfs)+' OUTPUT=raw.SNPs.vcf VALIDATION_STRINGENCY=SILENT'
f.write('wait\n'+mergeVCF+'\n')
f.close()
f2=open('vqsr.sh', 'w')

#-an DP for whole genome
vqsr1 = 'java  -Djava.io.tmpdir=/mnt/cluster2/tmp  -Xms1024m -Xmx4096m -jar  '+gatk+'GenomeAnalysisTK.jar -T VariantRecalibrator -R '+ref+' -input raw.SNPs.vcf '+ \
    ' -resource:hapmap,known=false,training=true,truth=true,prior=15.0 '+hapmap+ \
    ' -resource:omni,known=false,training=true,truth=true,prior=12.0 '+omni+ \
    ' -resource:1000G,known=false,training=true,truth=false,prior=10.0 '+onekg+ \
    ' -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 '+dbSNP+ \
    ' -an QD -an FS -an MQRankSum -an ReadPosRankSum '+\
    ' -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0'+ \
    ' -recalFile raw.SNPs.recal -tranchesFile raw.SNPs.tranches -rscriptFile recalibrate_SNP_plots.R'
vqsr2 = 'java -Djava.io.tmpdir=/mnt/cluster2/tmp  -Xms1024m -Xmx4096m -jar  '+gatk+'GenomeAnalysisTK.jar -T ApplyRecalibration -R '+ref+' -input raw.SNPs.vcf '+ \
   '-mode SNP -recalFile raw.SNPs.recal -tranchesFile raw.SNPs.tranches -o recal.SNPs.raw.indel.vcf -ts_filter_level 90.0' #change from 99.0
vqsr3='java -Djava.io.tmpdir=/mnt/cluster2/tmp  -Xms1024m -Xmx4096m -jar  '+gatk+'GenomeAnalysisTK.jar -T VariantRecalibrator -R '+ref+ \
   ' -input recal.SNPs.raw.indel.vcf'+ \
   ' -resource:mills,known=true,training=true,truth=true,prior=12.0 '+knownindel+ \
   ' -an MQRankSum -an FS -an ReadPosRankSum -mode INDEL '+ \
   ' -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 '+\
   ' --maxGaussians 4 ' + \
   ' -recalFile recalibrate_INDEL.recal -tranchesFile recalibrate_INDEL.tranches '+\
   ' -rscriptFile recalibrate_INDEL_plots.R '
vqsr4='java -Djava.io.tmpdir=/mnt/cluster2/tmp  -Xms1024m -Xmx4096m -jar  '+gatk+'GenomeAnalysisTK.jar -T ApplyRecalibration -R '+ref+ \
    ' -input recal.SNPs.raw.indel.vcf -mode INDEL --ts_filter_level 90.0 '\
    ' -recalFile recalibrate_INDEL.recal -tranchesFile recalibrate_INDEL.tranches '+\
    ' -o recalibrated_variants.vcf' 
phase1 ='java -Djava.io.tmpdir=/mnt/cluster2/tmp  -Xms1024m -Xmx4096m -jar  '+gatk+'GenomeAnalysisTK.jar -T PhaseByTransmission -R '+ref+\
    '-V recalibrated_variants.vcf'+ \
    '-ped'+fam +'-o phased1.vcf'
phase2='java -Djava.io.tmpdir=/mnt/cluster2/tmp  -Xms1024m -Xmx4096m -jar  '+gatk+'GenomeAnalysisTK.jar -T ReadBackedPhasing -R '+ref+' '+bams+\
     '-V phased1.vcf -o phased2.vcf -respectPhaseInput --phaseQualityThresh 20.0'
beagle1 ='java -Xmx4g -Djava.io.tmpdir=/mnt/cluster2/tmp  -Xms1024m -Xmx4096m -jar  '+gatk+'GenomeAnalysisTK.jar -L 20 -R '+ref+' -T ProduceBeagleInput '+\
      '-V phased2.vcf -o beagle.in'
#run beagle
#beagle output uncompress
beagle2='java -Djava.io.tmpdir=/mnt/cluster2/tmp -Xmx4000m  -Xms1024m -Xmx4096m -jar  '+gatk+'GenomeAnalysisTK.jar -R '+ref+' -T BeagleOutputToVCF '+\
      ' -V phased2.vcf '+\
      ' -beagleR2:BEAGLE /myrun.beagle_output.r2 '+\
      ' -beaglePhased:BEAGLE /myrun.beagle_output.phased '+\
      ' -beagleProbs:BEAGLE /myrun.beagle_output.gprobs '+\
      ' -o output_vcf.vcf'

# hardSNP1 = 'java -Djava.io.tmpdir=/mnt/cluster2/tmp -Xmx4000m  -Xms1024m -Xmx4096m -jar '+gatk+'GenomeAnalysisTK.jar '+\
    # '-T SelectVariants -R '+ref+' -V raw.SNPs.vcf -L 20 -selectType SNP -o raw_snpsonly.vcf '
	
# hardSNP2 = 'java -Djava.io.tmpdir=/mnt/cluster2/tmp -Xmx4000m  -Xms1024m -Xmx4096m -jar '+gatk+'GenomeAnalysisTK.jar '+\
    # '-T VariantFiltration -R '+ref+' -V raw_snpsonly.vcf '+ \
    # '--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0" '+\
    # '--filterName "hard_snp_filter" '+ \
    # '-o filtered_snps.vcf '

# hardIndel1 ='java -Djava.io.tmpdir=/mnt/cluster2/tmp -Xmx4000m  -Xms1024m -Xmx4096m -jar '+gatk+'GenomeAnalysisTK.jar '+\
    # '-T SelectVariants -R '+ref+' -V raw.SNPs.vcf -L 20 -selectType INDEL -o raw_indelsonly.vcf' 
# hardIndel2 ='java -Djava.io.tmpdir=/mnt/cluster2/tmp -Xmx4000m  -Xms1024m -Xmx4096m -jar '+gatk+'GenomeAnalysisTK.jar '+\
    # '-T VariantFiltration -R '+ref+' -V raw_indelsonly.vcf '+\
    # '--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"'+ \
    # '--filterName "hard_indel_filter" -o filtered_indels.vcf'

#wget -O dbSnp.vcf.gz ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz
#zcat dbSnp.vcf.gz > dbSnp.vcf
passOnly = '/mnt/cluster/tools/vcftools_0.1.12a/bin/vcftools --vcf recalibrated_variants.vcf --remove-filtered-all --recode --stdout > filtered.vcf'
#--bed SeqCap_EZ_Exome_v3_capture_fix.bed
# Annotate ID field using dbSnp
bedFilter ='java -Djava.io.tmpdir=/mnt/cluster2/tmp -Xmx4000m  -Xms1024m -Xmx4096m -jar '+gatk+'GenomeAnalysisTK.jar '+\
    '-T VariantFiltration -R '+ref+' --filterNotInMask --mask SeqCap_EZ_Exome_v3_capture_fix.bed -V filtered.vcf -o filtered.bed.vcf '   
passOnlyBed = '/mnt/cluster/tools/vcftools_0.1.12a/bin/vcftools --vcf filtered.bed.vcf --remove-filtered-all --recode --stdout > filtered.bed.mask.vcf'

#snpsift = 'java -Xmx4g  -Djava.io.tmpdir=/mnt/cluster2/tmp  -Xms1024m -Xmx4096m -jar  '+snpeff+'SnpSift.jar annotate -v /mnt/cluster2/gatk_resource2.8/dbSnp.vcf.gz recalibrated_variants.vcf > snpsift.vcf'
snpeff1 = 'java -Xmx4g  -Djava.io.tmpdir=/mnt/cluster2/tmp  -Xms1024m -Xmx4096m -jar  /mnt/san/cluster/tools/snpEff/snpEff.jar eff -c /mnt/san/cluster/tools/snpEff/snpEff.config -s filtered.summary.html -v GRCh37.74 filtered.vcf > snpEff.filtered.vcf'
snpeff2 = 'java -Xmx4g  -Djava.io.tmpdir=/mnt/cluster2/tmp  -Xms1024m -Xmx4096m -jar  /mnt/san/cluster/tools/snpEff/snpEff.jar eff -c /mnt/san/cluster/tools/snpEff/snpEff.config -s filteredBed.summary.html -v GRCh37.74 filtered.bed.mask.vcf > snpEff.filtered.bed.vcf'

# snpeff2 = 'java -Xmx4g -Djava.io.tmpdir=/mnt/cluster2/tmp  -Xms1024m -Xmx4096m -jar  '+gatk+'GenomeAnalysisTK.jar  -R '+ref+\
          # '-T VariantAnnotator -A SnpEff --variant snpEff_input.vcf --snpEffFile snpEff_output.vcf -o out.vcf\n'
#f2.write('\n'.join([allstats, call, vqsr1, vqsr2, vqsr3, vqsr4])+'\n')#
f2.write('\n'.join([vqsr1, vqsr2, vqsr3, vqsr4, passOnly, bedFilter, passOnlyBed, snpeff1, snpeff2, 'python ti_tv.py'])+'\n')
#f2.write('\n'.join([vqsr3, vqsr4, passOnly,bedFilter])+'\n')
#f2.write('\n'.join([ passOnlyBed, snpeff2])+'\n')
#f.write(allstats+'\n')
f2.close()

