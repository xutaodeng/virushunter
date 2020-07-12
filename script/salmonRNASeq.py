#!/usr/bin/env python
from collections import defaultdict
import os.path
import os
import string
import sys
from os import walk
from os import listdir
# step 1: install Salmon salmon-0.10.0_linux_x86_64, fastqc and AdapterRemoval
# step 2: download human reference GRCh38 CDS and ncRNA from ensembl
	#  https://uswest.ensembl.org/info/data/ftp/index.html
# step 3: run fastqc and adaptor removal
# step 4: process CDS sequence to remove first 5 codons 
# step 5: construct 2 reference sequence file: A. combine nc with CDS unaltered. B. combien nc with codon removed
# step 6: generate salmon index with refA and refB
	# Trim 15 codons
	# fa15codon.py Homo_sapiens.GRCh38.cds.all.fa Homo_sapiens.GRCh38.cds.no15Codons.fa
	# BSI\306307@bsidna1:/mnt/cluster/tools> cat Homo_sapiens.GRCh38.ncrna.fa Homo_sapiens.GRCh38.cds.no15Codons.fa > Homo_sapiens.GRCh38.nc.cds.no15Codons.fa
	# BSI\306307@bsidna1:/mnt/cluster/tools> cat Homo_sapiens.GRCh38.ncrna.fa Homo_sapiens.GRCh38.cds.all.fa > Homo_sapiens.GRCh38.nc.cds.all.fa
	#/mnt/cluster/xdeng/tools/salmon-0.10.0_linux_x86_64/bin/salmon index -t Homo_sapiens.GRCh38.nc.cds.all.fa  -i Homo_sapiens.GRCh38.cds.all.fa._index --type quasi -k 31

	#/mnt/cluster/xdeng/tools/salmon-0.10.0_linux_x86_64/bin/salmon index -t Homo_sapiens.GRCh38.nc.cds.no15Codons.fa  -i Homo_sapiens.GRCh38.cds.no15Codons.fa._index --type quasi -k 31
# step 7. quant fastq using Salmon with refA and refB separately
# step 8. Summarize TPM and ReadCount from two results, sumarize nc and coding sequences
	# final output file is codon.TPM and codon.reads with refB. full.TPM and full.reads with refA

servers=['bsidna3','bsidna4','bsidna5','bsidna6','bsidna7','bsidna8','bsidna9','bsidna10', 'bsidna11',\
'bsidna12','bsidna13','bsidna14','bsidna15','bsidna16','bsidna17','bsidna18','bsidna19', \
'bsidna21','bsidna22', 'bsidna23','bsidna25','bsidna26', 'bsidna27','bsidna28','bsidna29', 'bsidna30', \
'bsidna32', 'bsidna33', 'bsidna35', 'bsidna36']

nservers=len(servers)

# servers=['bsidna15','bsidna16','bsidna17','bsidna18','bsidna19','bsidna20',\
# 'bsidna21','bsidna22', 'bsidna27','bsidna28','bsidna29', 'bsidna30']


#servers=['bsidna32']


# /mnt/cluster/xdeng/tools/FastQC/fastqc
# /mnt/cluster/xdeng/tools/adapterremoval-2.1.7/build/AdapterRemoval



samples=['SRR2033082',\
'SRR2033083',\
'SRR2033084',\
'SRR2033085',\
'SRR2033086',\
'SRR2033087',\
'SRR2033088',\
'SRR2033089',\
'SRR2033090',\
'SRR2033091',\
'SRR2033092',\
'SRR2033093']

samples=['6h1-foot_S14_L002_R1_001.fastq', '6h2-foot_S16_L002_R1_001.fastq', '6h2-mRNA_S15_L002_R1_001.fastq', '6h1-mRNA_S13_L002_R1_001.fastq']

wd = os.path.abspath(os.path.dirname('.')).replace('san2', 'cluster2')

def trim2():
	f=open('trim.sh', 'w')
	f0=open('fastqc.sh', 'w')
	f1=open('ncRNA.sh', 'w')
	f2=open('salmon.sh', 'w')
	f3=open('salmon_fullLength.sh', 'w')
	f4=open('paste.sh', 'w')
	f5=open('summarize.sh', 'w')
	job=0
	p2,p3=[],[]
	f0.write('/mnt/cluster/xdeng/tools/FastQC/fastqc '+' '.join(samples)+'\n')
	for s in samples:
		serverTag = 'ssh '+servers[job%nservers]
		
		f.write(serverTag+' "cd '+wd+' && /mnt/cluster/xdeng/tools/adapterremoval-2.1.7/build/AdapterRemoval --file1 '+s+' --basename ' \
		+s+' --qualitymax 100 --trimns --trimqualities --gzip "& \n')

		f2.write(serverTag+' "cd '+wd+'&& /mnt/cluster/xdeng/tools/salmon-0.10.0_linux_x86_64/bin/salmon quant -i /mnt/cluster/tools/Homo_sapiens.GRCh38.cds.no15Codons.fa._index -l A -r '+ \
			s+' -o '+s+'_out '+'"& \n')
		
		f3.write(serverTag+' "cd '+wd+'&& /mnt/cluster/xdeng/tools/salmon-0.10.0_linux_x86_64/bin/salmon quant -i /mnt/cluster/tools/Homo_sapiens.GRCh38.cds.all.fa._index -l A -r '+ \
			s+' -o '+s+'_full_out '+'"& \n')
		p2.append(s+'_out/quant.sf')
		p3.append(s+'_full_out/quant.sf')
		job+=1
		if job%nservers==0: 
			f.write('wait\n')
			f2.write('wait\n')
			f3.write('wait\n')
	f4.write('paste '+' '.join(p2)+' > '+'codon.sf \n')
	f4.write('paste '+' '.join(p3)+' > '+'full.sf \n')
	f5.write('salmon_summarize.py\n')
	f.close()
	f2.close()
	f3.close()
	f4.close()
	f5.close()
	f0.close()
	
def trim():
	f=open('trim.sh', 'w')
	f0=open('fastqc.sh', 'w')
	f1=open('ncRNA.sh', 'w')
	f2=open('salmon.sh', 'w')
	f3=open('salmon_fullLength.sh', 'w')
	f4=open('paste.sh', 'w')
	f5=open('summarize.sh', 'w')
	job=0
	p2,p3=[],[]
	f0.write('/mnt/cluster/xdeng/tools/FastQC/fastqc '+' '.join(samples)+'\n')
	for s in samples:
		serverTag = 'ssh '+servers[job%nservers]
		
		f.write(serverTag+' "cd '+wd+' && /mnt/cluster/xdeng/tools/adapterremoval-2.1.7/build/AdapterRemoval --file1 '+s+'/'+s+'.fastq --basename ' \
		+s+' --trimns --trimqualities --gzip "& \n')

		f2.write(serverTag+' "cd '+wd+'&& /mnt/cluster/xdeng/tools/salmon-0.10.0_linux_x86_64/bin/salmon quant -i /mnt/cluster/tools/Homo_sapiens.GRCh38.cds.no15Codons.fa._index -l A -r '+ \
			'./fastq/'+s+'.truncated.gz -o '+s+'_out '+'"& \n')
		
		f3.write(serverTag+' "cd '+wd+'&& /mnt/cluster/xdeng/tools/salmon-0.10.0_linux_x86_64/bin/salmon quant -i /mnt/cluster/tools/Homo_sapiens.GRCh38.cds.all.fa._index -l A -r '+ \
			'./fastq/'+s+'.truncated.gz -o '+s+'_full_out '+'"& \n')
		p2.append(s+'_out/quant.sf')
		p3.append(s+'_full_out/quant.sf')
		job+=1
		if job%nservers==0: 
			f.write('wait\n')
			f2.write('wait\n')
			f3.write('wait\n')
	f4.write('paste '+' '.join(p2)+' > '+'codon.sf \n')
	f4.write('paste '+' '.join(p3)+' > '+'full.sf \n')
	f5.write('salmon_summarize.py\n')
	f.close()
	f2.close()
	f3.close()
	f4.close()
	f5.close()
	f0.close()

trim2()

sys.exit()


#

def renamefastq():
	f=open('sampleInfo.txt', 'r')
	of=open('renamefastq.sh', 'w')
	for line in f:
		infile, outfile=line.strip().split()[0:2]
		infile=infile+'/s_1_1_sequence.txt'
		outfile='./fastq/'+outfile+'.fastq'
		of.write('mv '+infile+' '+outfile+'\n')
	f.close()
	of.close()

def readGroup():
	Group=defaultdict(list)
	try:
		f=open('samples_group.txt', 'r')
		for line in f:
			sample, group = line.strip().split()
			Group[group].append(sample)
	except:
		for key in seeds.keys():
			Group[key].append(key)
	print 'groups', Group
	return Group

def readContrast():
	Contrast=[]
	try:
		f=open('contrast.txt', 'r')
		for line in f:
			C1, C2= line.strip().split()
			Contrast.append([C1, C2])
	except:
		pass
	print 'contrast', Contrast
	return Contrast


def readdir(seedfile='./fastq/samples.txt'):
	seeds={}
	# ls -1d */
	f=open(seedfile, 'r')
	for line in f:
		directory=line.strip().rstrip('/')
		#print directory
		files=[wd+'/fastq/'+directory+'/'+x for x in listdir('./fastq/'+directory) if x.endswith('.gz')]
	# for (dirpath, dirnames, filenames) in walk(mypath):
		# files=[wd+'/'+dirpath[2:]+'/'+x for x in filenames]
		seeds[directory]=files
	#print seeds
	f.close()
	# print seeds.keys()
	# print '--------------------------------------------------'
	# for seed in seeds.keys():
		# print seed, seeds[seed]
	return seeds
		
	# seeds=defaultdict(list)
	# for line in f:
		# key2=line.strip().rsplit('.', 1)[0]
		# key1=line.strip().rsplit('_',4)[0]#split('_')[0]
		# if len(key1)<=len(key2):
			# key=key1
		# else:
			# key=key2
		# ss=line.strip()
		# if line.strip().endswith('fasta'): ss = ss.replace('fasta', 'fastq')
		# elif line.strip().endswith('sam'): ss = ss.replace('sam', 'fastq')
		# elif line.strip().endswith('bam'): ss = ss.replace('bam', 'fastq')
		# seeds[key].append(wd+'/fastq/'+ss)
	# f.close()
	# return seeds

def readdir_Pair(seedfile='fastq/samples.txt'):
	seeds={}
	# ls -1d */
	f=open(seedfile, 'r')
	for line in f:
		directory=line.strip()
		#print directory
		files=[wd+'/fastq/'+directory+'/'+x for x in listdir('./fastq/'+directory) if x.endswith('.gz')]
		files1=[x for x in files if '_R1_' in x]
		files2=[x for x in files if '_R2_' in x]
		files1.sort()
		files2.sort()
		print files1, files2
		print len(files1), len(files2)
		assert len(files1)==len(files2)
		for i in xrange(len(files1)):
			assert files1[i].replace('_R1_', '_R2_')==files2[i]
		seeds[directory]=(files1, files2)

	f.close()
	#print seeds.keys()
	print '--------------------------------------------------'
	for seed in seeds.keys():
		print seed, seeds[seed]
	return seeds
		
def readSeeds2(seedfile='fastq/samples.txt'):
	seeds=defaultdict(list)
	f=open(seedfile, 'r')
	for line in f:
		key2=line.strip().rsplit('.', 1)[0]
		key1=line.strip().rsplit('_',1)[0]#split('_')[0]
		#print line
		# key1=line.strip().split('_',1)[0]#split('_')[0]
		# if len(key1)<=len(key2):
			# key=key1
		# else:
			# key=key2
		key=key1
		ss=line.strip()
		if line.strip().endswith('fasta'): ss = ss.replace('fasta', 'fastq')
		elif line.strip().endswith('sam'): ss = ss.replace('sam', 'fastq')
		elif line.strip().endswith('bam'): ss = ss.replace('bam', 'fastq')
		seeds[key].append(wd+'/fastq/'+ss)
	f.close()
	return seeds

def readSeeds3(seedfile='fastq/samples.txt'):
	seeds=defaultdict(list)
	f=open(seedfile, 'r')
	for line in f:
		key='_'.join(line.strip().rsplit('_')[0:2])#split('_')[0]
		ss=line.strip()
		if line.strip().endswith('fasta'): ss = ss.replace('fasta', 'fastq')
		elif line.strip().endswith('sam'): ss = ss.replace('sam', 'fastq')
		elif line.strip().endswith('bam'): ss = ss.replace('bam', 'fastq')
		seeds[key].append(wd+'/fastq/'+ss)
	f.close()
	return seeds

def readSeeds4(seedfile='fastq/samples.txt'):
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

	
def dedup(pair):
	f = open('dedup.sh', 'w')
	job=0
	for (key, fqfiles) in seeds.items():
		if pair: 
			#print fqfiles
			fq1, fq2=fqfiles
			fqs=fq1+fq2
		else:
			fqs=fqfiles
		for fq in fqs:
			serverTag = 'ssh '+servers[job%nservers]
			f.write(serverTag + ' '+dirscr+'dedup.py '+fq+' '+fq+'.dup &\n')
			job+=1
			if job%nservers==0:  f.write('wait\n')
			
def tophat(pair):
	f = open('tophat.sh', 'w')
	job=0
	for (key, fqfiles) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		
		if not pair: f.write(serverTag+' "'+tophatPath+' -o '+wd+'/'+key+'_t_out '+genome+' '.join(fqfiles)+' "&\n')
		if pair: 
			fq1, fq2=fqfiles
			print 'haha', fq1, fq2
			if type(fq1) is list:
				f.write(serverTag+' "'+tophatPath+' -o '+wd+'/'+key+'_t_out '+genome+','.join(fq1)+' '+','.join(fq2)+' "&\n')
			else:
				f.write(serverTag+' "'+tophatPath+' -o '+wd+'/'+key+'_t_out '+genome+' '+fq1+' '+fq2+' "&\n')
		job+=1
		if job%nservers==0:  f.write('wait\n')
	f.close()

def samsort():
	f = open('samsort.sh', 'w')
	job=0
	i=0
	for (key, fqfiles) in seeds.items():
		base=wd+'/'+key+'_t_out/accepted_hits'
		infile = base+'_'+key+'.sam'
		outfile =infile+'.sort'
		serverTag = 'ssh '+servers[job%nservers]
		f.write(serverTag+' /mnt/cluster/xdeng/tools/samtools-1.2/samtools sort '+infile+' '+outfile+' &\n')
		job+=1
		if job%nservers==0: f.write('wait\n')
		i+=1
	f.write('wait\n')
	f.close()

# def mpileup():
	# bams=[]
	# f = open('pileup.sh', 'w')
	# f.write('/mnt/cluster/xdeng/tools/samtools-1.2/samtools mpileup  -d 1000000 -C50 -Df '+genomefa+' ')
	# job=0
	# for (key, fqfiles) in seeds.items():
		# serverTag = 'ssh '+servers[job%nservers]
		# base=wd+'/'+key+'_t_out/accepted_hits'
		# bams.append(base+'.sort.bam')
		# serverTag = 'ssh '+servers[job%nservers]
		# job+=1

	# f.write(' '.join(bams)+' > all.pileup\n')
	# f.close()


def mpileup_precise():
	bams=[]
	f = open('pileup_precise.sh', 'w')
	f.write('/mnt/cluster/xdeng/tools/samtools-1.2/samtools mpileup -BQ0 -d10000000 -f '+genomefa+' ')
	for (key, fqfiles) in seeds.items():
		base=wd+'/'+key+'_t_out/accepted_hits'
		bams.append(base+'_'+key+'.sam.sort.bam')
	f.write(' '.join(bams)+' > all.precise.pileup\n')
	f.close()

def finalpileup():
	bams=[]
	f = open('finalpileup.sh', 'w')
	#f.write(dirscr+'pileup.py all.pileup pileup.sh count.pileup percent.pileup '+ref+'\n')
	f.write(dirscr+'pileup.py all.precise.pileup pileup_precise.sh count.pileup percent.pileup '+genomefa+'\n')
	f.close()

def readchrosize(organism):
    directory=os.path.dirname(sys.argv[0])
    chrofile=directory+'/ChromSizes/'+organism+'.chrom.sizes'
    print 'chrosize file', chrofile
    f=open(chrofile, 'r')
    chrosize={}
    for line in f:
        chro, size = line.strip().split()
        chro=chro.upper()
        if chro.startswith('CHRO'):
            chro=chro.replace('CHRO', '')
        elif chro.startswith('CHR'):
            chro=chro.replace('CHR', '')
        elif chro.startswith('CH'):
            chro=chro.replace('CH', '')
        size = int(size)
        chrosize[chro]=size
    f.close()
    print 'chro size', chrosize
    return chrosize, chrofile

def bam2bigwig(region):
	f0 = open('samindex.sh', 'w')
	f = open('bam2sam.sh', 'w')
	f2 = open('sam2wig.sh', 'w')
	#f3 = open('wig2bigwig.sh', 'w')
	job=0
	#if region is '', whole genome will be processed
	for (key, fqfiles) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		base=wd+'/'+key+'_t_out/accepted_hits'
		f0.write(serverTag+' "/mnt/cluster/xdeng/tools/samtools-1.3/samtools index '+base+'.bam" &\n') 
		f.write(serverTag+' "/mnt/cluster/xdeng/tools/samtools-1.3/samtools view -h -o '+base+'_'+key+'.sam '+base+'.bam "'+region+'" "&\n')
		#f.write(serverTag+' "/mnt/cluster/xdeng/tools/samtools-1.3/samtools view -h -o '+base+'.sam '+base+'.bam "&\n')
		f2.write(serverTag+' "/mnt/cluster/xdeng/script/sam2wig.py -n 0 -o '+chrosize+' '+base+'_'+key+'.sam "&\n') #no normalize
		# f3.write(serverTag+' "/mnt/san/cluster/bsidna1_local/ChIPseeqer-2.1/dist/SCRIPTS/wigToBigWig -clip '+ \
				# base+'.sam.wig '+chrofile+' '+base+'.sam.wig.bw"&\n')
		job+=1
		if job%nservers==0:  
			f.write('wait\n')
			f2.write('wait\n')
			f0.write('wait\n')
	f.close()
	f2.close()
	f0.close()
	#f3.close()
	#samtools view -h -o out.sam in.bam
	
#http://uswest.ensembl.org/info/data/ftp/index.html
#/mnt/cluster/mm10/Mus_musculus/UCSC/mm10/Sequence/AbundantSequences/nc
def bowtie_NC(pair):
	try: os.mkdir(wd+'/NC/')
	except: pass
	f = open('bowtie_NC.sh', 'w')
	f2 = open('bowtie_NC_count.sh', 'w')
	job=0
	for (key, fqfiles) in seeds.items():
		sam=wd+'/NC/'+key+'.sam'
		samC=wd+'/NC/'+key+'.sam.count'
		serverTag = 'ssh '+servers[job%nservers]
		#print 'fqfiles', fqfiles
		
		if not pair: 
			fqs=[' -U ' +fq for fq in fqfiles]
			f.write(serverTag+' "'+bowtie2Path+'  --quiet --local --no-hd -p 7 -x '+nc+' '+' '.join(fqs)+' -S '+sam+' &\n')
		if pair: 
			fq1, fq2=fqfiles
			f.write(serverTag+' "'+bowtie2Path+'  --quiet --local --no-hd -p 7 -x '+nc+' -1 '+fq1+' -2 '+fq2+' -S '+sam+' &\n')
		f2.write(serverTag+' '+dirscr+'/non_codingRNA.py '+sam+' '+ncfa+' '+samC+' &\n')
		job+=1
		if job%nservers==0:  
			f.write('wait\n')
			f2.write('wait\n')
	f.close()
	f2.close()

def cufflinks():
	f = open('cufflinks.sh', 'w')
	job=0
	for (key, fqfiles) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		f.write(serverTag+' '+cufflinksPath+' -o '+wd+'/'+key+'_c_out '+wd+'/'+key+'_t_out/accepted_hits.bam &\n')
		job+=1
		if job%nservers==0:  f.write('wait\n')
	f.close()
def readCount():
	of=open('read.count', 'w')
	of.close()
	f = open('readCount.sh', 'w')
	job=0
	for (key, fqfiles) in seeds.items():
		serverTag = ' ' #'ssh '+servers[job%nservers]
		f.write(dirscr+'/readCount.py '+wd+'/'+key+'_t_out/accepted_hits.bam '+wd+'/'+key+' \n')
		#job+=1
		#if job%nservers==0:  f.write('wait\n')
	f.close()
	
def cuffmerge():
	f = open('cuffmerge.sh', 'w')
	f2 = open('assemblies.txt', 'w')
	for (key, fqfiles) in seeds.items():
		f2.write('./'+key+'_c_out/transcripts.gtf\n')
	f2.close()
	f.write(cuffmergePath+'\n')
	f.close()

def cuffquant():
	f = open('cuffquant.sh', 'w')
	job=0
	for (key, fqfiles) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		f.write(serverTag+' '+cuffquantPath+' -o '+wd+'/'+key+'_q_out '+' '+wd+'/merged_asm/merged.gtf '+wd+'/'+key+'_t_out/accepted_hits.bam &\n')
		job+=1
		if job%nservers==0:  f.write('wait\n')
	f.close()

def cuffnorm():
	f = open('cuffnorm.sh', 'w')
	f.write(cuffnormPath+' '+wd+'/merged_asm/merged.gtf ')
	for (key, fqfiles) in seeds.items():
		f.write(wd+'/'+key+'_q_out/abundances.cxb ')
	f.write('\n')
	f.close()

def cuffdiff(Group, Contrast): #two groups set for comparison and output label
	f = open('cuffdiff.sh', 'w')
	f2=open('cuffdiff_combine_result.sh', 'w')
	job=0
	labels=[]
	for con in Contrast:
		C1key, C2key=con
		C1, C2=Group[C1key], Group[C2key]
		label=C1key+'_'+C2key
		serverTag = 'ssh '+servers[job%nservers]
		cmd=serverTag+' '+cuffdiffPath+' -b '+genomefa+' -L C1,C2 -u '+wd+'/merged_asm/merged.gtf -o '+wd+'/'+label+'/diff_d_out '
		labels.append(label)
		g1,g2=[],[]
		for (key, fqfiles) in seeds.items():
			if key in C1: 
				g1.append(wd+'/'+key+'_t_out/accepted_hits.bam')
			if key in C2: 
				g2.append(wd+'/'+key+'_t_out/accepted_hits.bam')
		cmd+=','.join(g1)+' '+','.join(g2)
		f.write(cmd+' & \n')
		job+=1
		if job%nservers==0:  f.write('wait\n')
	f.close()
	f2.write('RNASeq_combine.py '+' '.join(labels)+'\n')
	f2.close()


if __name__ == "__main__":
	organism = sys.argv[1]
	pair=False
	wd = os.path.abspath(os.path.dirname('.'))
	dirscr='/mnt/cluster/xdeng/script/'
	if not wd.startswith('/mnt/'): wd='/mnt'+wd
	wd=wd.replace('/san2/', '/cluster2/')
	
	if organism == 'mouse':
		genome = "/mnt/cluster/mm10/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome "
		gtf=" /mnt/cluster/mm10/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf "
		genomefa=" /mnt/cluster/mm10/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome.fa "
		nc = "/mnt/cluster/mm10/Mus_musculus/UCSC/mm10/Sequence/AbundantSequences/nc "
		ncfa = "/mnt/cluster/mm10/Mus_musculus/UCSC/mm10/Sequence/AbundantSequences/nc.fasta "
		chrosize = " mm10 "
		chrofile='mnt/cluster/xdeng/script/ChromSizes/mm10.chrom.sizes'
	else:
		genomefa=" /mnt/cluster/hg19/Sequence/Bowtie2Index/genome.fa "
		gtf ="/mnt/cluster/hg19/Annotation/Genes/genes.gtf "
		genome ="/mnt/cluster/hg19/Sequence/Bowtie2Index/genome "
		nc ="/mnt/cluster/hg19/Sequence/AbundantSequences/nc "
		ncfa ="/mnt/cluster/hg19/Sequence/AbundantSequences/nc.fasta "
		chrosize = " hg19 "
		chrofile='mnt/cluster/xdeng/script/ChromSizes/hg19.chrom.sizes'

	#strandSpecific=" --library-type fr-firststrand "
	strandSpecific=" "
	bowtie2Path='/mnt/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2 '
	tophatPath="/mnt/cluster/tools/tophat-2.0.14.Linux_x86_64/tophat -p 8 --no-coverage-search -G "+gtf +strandSpecific# reads_1.fq reads_2.fq
	cufflinksPath="/mnt/cluster/tools/cufflinks-2.2.1.Linux_x86_64/cufflinks -p 8"
	cuffmergePath="/mnt/cluster/tools/cufflinks-2.2.1.Linux_x86_64/cuffmerge -p 8 -g "+gtf+" -s "+genomefa+" assemblies.txt"
	cuffquantPath="/mnt/cluster/tools/cufflinks-2.2.1.Linux_x86_64/cuffquant -p 8 "
	cuffnormPath="/mnt/cluster/tools/cufflinks-2.2.1.Linux_x86_64/cuffnorm -p 8 "
	cuffdiffPath="/mnt/cluster/tools/cufflinks-2.2.1.Linux_x86_64/cuffdiff -p 8"
	#renamefastq()
	#sys.exit()
	#seeds=readSeeds4()
	#seeds=readSeeds4()
	seeds=readSeeds2()
	for key in seeds.keys():
		print key, seeds[key]

	dedup(pair)
	tophat(pair)
	readCount()
	cufflinks()
	cuffmerge()
	cuffquant()
	cuffnorm()
	bowtie_NC(pair)
	samsort()
	finalpileup()
	mpileup_precise()
	
	
	# Group= readGroup()
	# Contrast= readContrast()
	# cuffdiff(Group, Contrast)
	# bam2bigwig('chr3:22076652-22216594')
	
	sf=open('pipeline_run.sh', 'w')
	sf.write('source tophat.sh \nwait\n')
	sf.write('source readCount.sh \nwait\n')
	sf.write('source bowtie_NC.sh \nwait\n')
	sf.write('source bowtie_NC_count.sh \nwait\n')
	sf.write('source cufflinks.sh \nwait\n')
	sf.write('source cuffmerge.sh \nwait\n')
	sf.write('source cuffquant.sh \nwait\n')
	sf.write('source cuffnorm.sh \nwait\n')
	sf.write('source cuffdiff.sh \n')
	sf.write('source cuffdiff_combine_result.sh \n')
	sf.write('source samsort.sh\nwait\n')
	sf.write('source pileup_precise.sh\nwait\n')
	sf.write('source finalpileup.sh\nwait\n')
	sf.close()
	print 'organism', organism
	print 'pair', pair
	print 'strand specific=',strandSpecific
	print 'wd', wd