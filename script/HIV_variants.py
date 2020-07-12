#!/usr/bin/env python
from collections import defaultdict
import os.path
import os
import string
import sys

servers=['bsidna2','bsidna3', 'bsidna4','bsidna5','bsidna6','bsidna7','bsidna8','bsidna9','bsidna10','bsidna11',\
'bsidna12','bsidna13','bsidna14','bsidna15','bsidna16','bsidna17','bsidna18','bsidna19','bsidna20',\
'bsidna21','bsidna22', 'bsidna23', 'bsidna24','bsidna25', 'bsidna26','bsidna27','bsidna28','bsidna29', 'bsidna30', \
'bsidna32','bsidna33', 'bsidna34', 'bsidna35','bsidna36']

# Non-Nextera

# 1) Assembling paired ends
# 2) Quality filtration (based on scores in FastQ): end trimming / deletion / low complexity
# 3) Deduplication -> List of unique variants with relative frequencies*
# 4) No dedup!
# 5) "Error" filtration (<1%)
# 6) Between-sample distance using reference
# 7) Nucleotide variant analyses from both data sets (w/ respect to HXB2 reference)
# 8) Take all protein sequences (from both data sets) -> only "viable" sequences with no premature stop codons.
# 9) Protein variant analyses from both data sets (w/ respect to HXB2 reference)

# Nextera
# 1) Assembling paired ends
# 2) Quality filtration (based on scores in FastQ): end trimming / deletion / low complexity
# 3) Nucleotide variant analysis (w/ respect to HXB2 reference)
# 4) Protein variant analysis (w respect to HXB2 reference)

dirscr='/mnt/cluster/xdeng/script/'
bowtieHumanIndex=' /mnt/cluster/hg19/Sequence/Bowtie2Index/genome '
nservers=len(servers)

def bowtieBuild(ref, bowindex):
	os.system('/mnt/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build '+ref+' '+ bowindex)

def readSeeds2(seedfile='fastq/samples.txt'):
	seeds=defaultdict(list)
	f=open(seedfile, 'r')
	for line in f:
		key2=line.strip().rsplit('.', 1)[0]
		key1=line.strip().split('_')[0]
		if len(key1)<len(key2):
			key=key1
		else:
			key=key2
		ss=line.strip()
		if line.strip().endswith('fasta'): ss = ss.replace('fasta', 'fastq')
		elif line.strip().endswith('sam'): ss = ss.replace('sam', 'fastq')
		elif line.strip().endswith('bam'): ss = ss.replace('bam', 'fastq')
		seeds[key].append(ss)
	f.close()
	return seeds

def readSeedsNew(seedfile='fastq/samples.txt'):
	seeds=defaultdict(list)
	seeds2=defaultdict(set)
	f=open(seedfile, 'r')
	for line in f:
		key2=line.strip().rsplit('.', 1)[0]
		key1=line.strip().split('_')[0]
		if len(key1)<len(key2):
			key=key1
		else:
			key=key2
		parts = key.split('-')
		newkey=parts[0]+'-'+parts[1][:-1]+parts[2]
		seeds2[newkey].add(key)
		ss=line.strip()
		if line.strip().endswith('fasta'): ss = ss.replace('fasta', 'fastq')
		elif line.strip().endswith('sam'): ss = ss.replace('sam', 'fastq')
		elif line.strip().endswith('bam'): ss = ss.replace('bam', 'fastq')
		seeds[key].append(ss)
	f.close()
	return seeds, seeds2


def catFq():
	f = open('catFq.sh', 'w')
	job=0
	seeds3=defaultdict(list)
	for (newkey, keys) in seeds2.items():
		fq1, fq2=[],[]
		for key in keys:
			fqfile1=seeds[key][0]
			fqfile2=seeds[key][1]
			fq1.append(wd+'/fastq/'+fqfile1)
			fq2.append(wd+'/fastq/'+fqfile2)
			if fqfile1[0].endswith('.gz'):
				cmd = ' zcat '
			else:
				cmd = ' cat '
			serverTag = 'ssh '+servers[job%nservers]
		f.write(serverTag+cmd+' '+' '.join(fq1)+' > '+wd+'/fastq/'+newkey+'_1.fq &\n')
		f.write(serverTag+cmd+' '+' '.join(fq2)+' > '+wd+'/fastq/'+newkey+'_2.fq &\n')
		seeds3[newkey]=[newkey+'_1.fq', newkey+'_2.fq']
		job+=1
	f.close()
	print len(seeds), seeds
	return seeds3


def dna2prot():
	f=open('dna2prot.sh', 'w')
	job=0
	for (key, fqfiles) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		f.write(serverTag+' '+dirscr+'dna2prot.py '+ wd+'/'+key+'.fa '+wd+'/'+key+'_prot.fa True &\n')
		f.write(serverTag+' '+dirscr+'dna2prot.py '+ wd+'/dedupDNA/'+key+'_dedup.fa '+wd+'/dedupProt/'+key+'_dedup_prot.fa True &\n')
		job+=1
	f.write('wait\n')
	f.close()

def mergePairFq():
	f = open('mergePairFq.sh', 'w')
	job=0
	for (key, fqfiles) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		fq = [wd+'/fastq/'+fqfile for fqfile in fqfiles]
		f.write(serverTag+' /mnt/cluster/xdeng/tools/FLASH-1.2.11/flash -o '+key+' -d '+wd+' '+' '.join(fq)+' &\n')
		job+=1
	f.write('wait\n')
	job=0
	for (key, fqfiles) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		f.write(serverTag+' cat '+wd+'/'+key+'.notCombined_1.fastq '+ wd+'/'+key+'.notCombined_2.fastq '+wd+'/'+key+'.extendedFrags.fastq > '+wd+'/'+key+'.fq &\n')
		job+=1
	f.close()
	
def mergeRunFq():
	f = open('mergeRunFq.sh', 'w')
	job=0
	for (key, fqfiles) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		fq = [wd+'/fastq/'+fqfile for fqfile in fqfiles]
		
		if fqfiles[0].endswith('.gz'):
			f.write(serverTag+' zcat '+' '.join(fq)+' > '+wd+'/'+key+'.fq &\n')
		else:
			f.write(serverTag+' cat '+' '.join(fq)+' > '+wd+'/'+key+'.fq &\n')
		job+=1
	f.write('wait\n')
	f.close()
	

def trimFq():
	f = open('trim2ends.sh', 'w')
	job=0
	for (key, fqfiles) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		f.write(serverTag+' '+dirscr+'trim2ends.py '+wd+'/'+key+'.fq '+wd+'/'+key+'.trim &\n')
		job+=1
	f.write('wait\n')
	f.close()

def fq2fa():
	f = open('fq2fa.sh', 'w')
	job=0
	for (key, fqfiles) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		f.write(serverTag+' '+dirscr+'fq2fa.py '+wd+'/'+key+'.trim '+wd+'/'+key+'.fa &\n')
		job+=1
	f.write('wait\n')
	f.close()

def dedup():
	f = open('dedup.sh', 'w')
	job=0
	try: os.mkdir(wd+'/dedupDNA/')
	except: pass
	try: os.mkdir(wd+'/dedupProt/')
	except: pass
	for (key, fqfiles) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		f.write(serverTag+' '+dirscr+'HIV_dedup.py '+wd+'/'+key+'.fa '+wd+'/dedupDNA/'+key+'_dedup.fa 0.0001 &\n')
		job+=1
	f.write('wait\n')
	f.close()


def mpileup():
	bams=[]
	f = open('pileup.sh', 'w')
	f2 = open('consensus.sh', 'w')
	f.write('/mnt/cluster/xdeng/tools/samtools-1.2/samtools mpileup  -d 1000000 -C50 -Df '+ref+' ')
	job=0
	for (key, fqfiles) in seeds.items():
		bams.append(key+'.sort.bam')
		serverTag = 'ssh '+servers[job%nservers]
		# f2.write(serverTag+' '+'/mnt/cluster/xdeng/tools/samtools-1.3/samtools mpileup -uf '+ref+' '+wd+'/'+ key+'.sort.bam  | ' + \
		       # ' /mnt/cluster/xdeng/tools/bcftools-1.3/bcftools view -cg - | /mnt/cluster/xdeng/tools/bcftools-1.3/vcfutils.pl vcf2fq > ' +\
			   # wd+'/'+key+'consensus.fq &\n')
		job+=1
		if job%nservers==0: f2.write('wait\n')
	f.write(' '.join(bams)+' > all.pileup\n')
	f.close()
	f2.close()

def mpileup_precise():
	bams=[]
	f = open('pileup_precise.sh', 'w')
	f.write('/mnt/cluster/xdeng/tools/samtools-1.2/samtools mpileup -BQ0 -d10000000 -f '+ref+' ')
	for (key, fqfiles) in seeds.items():
		bams.append(key+'.sort.bam')
	f.write(' '.join(bams)+' > all.precise.pileup\n')
	f.close()
	
def sam2wig():
	bams=[]
	f = open('sam2wig.sh', 'w')
	job=0
	for (key, fqfiles) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		f.write(serverTag+' '+dirscr+'sam2wig.py  '+wd+'/'+key+'.sam &\n')
		job+=1
		if job%nservers==0: f.write('wait\n')
	f.close()

def wig2svg():
	bams=[]
	f = open('wig2svg.sh', 'w')
	job=0
	for (key, fqfiles) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		f.write(serverTag+' '+dirscr+'wig2svg.py  '+wd+'/'+key+'.sam.wig '+wd+'/'+key+'.html &\n')
		job+=1
		if job%nservers==0: f.write('wait\n')
	f.close()

def sam2wig_multi():
	bams=[]
	f = open('sam2wig_multi.sh', 'w')
	job=0
	for (key, fqfiles) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		f.write(serverTag+' '+dirscr+'sam2wig.py  '+wd+'/'+key+'.multi.sam '+wd+'/'+key+'.multi.wig &\n')
		job+=1
		if job%nservers==0: f.write('wait\n')
	f.close()
	
def sam2wig_human():
	bams=[]
	f = open('sam2wig_human.sh', 'w')
	job=0
	for (key, fqfiles) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		f.write(serverTag+' '+dirscr+'sam2wig.py  '+wd+'/'+key+'.human.sam '+wd+'/'+key+'.human.wig &\n')
		job+=1
		if job%nservers==0: f.write('wait\n')
	f.close()

def finalpileup():
	bams=[]
	f = open('finalpileup.sh', 'w')
	#f.write(dirscr+'pileup.py all.pileup pileup.sh count.pileup percent.pileup '+ref+'\n')
	f.write(dirscr+'pileup.py all.precise.pileup pileup_precise.sh count.pileup percent.pileup '+ref+'\n')
	f.close()
	
	
def bowtie():
	#/mnt/san/cluster/tools/bowtie2-2.1.0/bowtie2-build ENV_new.fa ENV_new
	f = open('bowtie.sh', 'w')
	job=0
	i=0
	for (key, fqfiles) in seeds.items():
		infile = wd+'/'+key+'.trim'
		outfile = wd+'/'+key+'.sam'
		serverTag = 'ssh '+servers[job%nservers]
		f.write(serverTag+' /mnt/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2 -D 20 -R 3 -N 1 -L 10 -i S,1,0.50 --quiet --local --reorder -p 7 -x '+bowtieindexpath+' -U '+infile+' -S '+outfile+' &\n')
		job+=1
		if job%nservers==0: f.write('wait\n')
		i+=1
	f.write('wait\n')
	f.close()

def bowtie_multi():
	#/mnt/san/cluster/tools/bowtie2-2.1.0/bowtie2-build ENV_new.fa ENV_new
	f = open('bowtie_multi.sh', 'w')
	job=0
	i=0
	for (key, fqfiles) in seeds.items():
		infile = wd+'/'+key+'.trim'
		outfile = wd+'/'+key+'.multi.sam'
		serverTag = 'ssh '+servers[job%nservers]
		f.write(serverTag+' /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2 -D 20 -R 3 -N 1 -L 10 -i S,1,0.50 --quiet --local --reorder -p 7 -x '+bowtieindexpath2+' -U '+infile+' -S '+outfile+' &\n')
		job+=1
		if job%nservers==0: f.write('wait\n')
		i+=1
	f.write('wait\n')
	f.close()

def bowtie_human():
	#/mnt/san/cluster/tools/bowtie2-2.1.0/bowtie2-build ENV_new.fa ENV_new
	f = open('bowtie_human.sh', 'w')
	job=0
	i=0
	
	for (key, fqfiles) in seeds.items():
		infile = wd+'/'+key+'.trim'
		outfile = wd+'/'+key+'.human.sam'
		serverTag = 'ssh '+servers[job%nservers]
		f.write(serverTag+' /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2 --quiet --local --reorder -p 7 -x '+bowtieHumanIndex+' -U '+infile+' -S '+outfile+' &\n')
		job+=1
		if job%nservers==0: f.write('wait\n')
		i+=1
	f.write('wait\n')
	f.close()

def samsort():
	f = open('samsort.sh', 'w')
	job=0
	i=0
	for (key, fqfiles) in seeds.items():
		infile = wd+'/'+key+'.sam'
		outfile = wd+'/'+key+'.sort'
		serverTag = 'ssh '+servers[job%nservers]
		f.write(serverTag+' /mnt/cluster/xdeng/tools/samtools-1.2/samtools sort '+infile+' '+outfile+' &\n')
		job+=1
		if job%nservers==0: f.write('wait\n')
		i+=1
	f.write('wait\n')
	f.close()

if __name__ == "__main__":
	ref = sys.argv[1]
	bowindex=ref.split('.')[0]
	bowtieBuild(ref, bowindex)
	#bowtieBuild('References_Sequences.fasta', 'References_Sequences')
	wd = os.path.abspath(os.path.dirname('.'))
	ref=wd+'/'+ref
	bowtieindexpath=wd+'/'+bowindex
	#bowtieindexpath2=wd+'/References_Sequences'
	#genseedfile()
	#if pair: 
	seeds=readSeeds2()
	# seeds, seeds2=readSeedsNew()
	# print 'nsample', len(seeds2)
	# for newkey in seeds2.keys():
		# print newkey, seeds2[newkey]
	# else:
		# seeds=readSeeds1()
	if not wd.startswith('/mnt/'): wd='/mnt'+wd
	for key in seeds.keys():
		print key, seeds[key]
	print 'num samples:', len(seeds)

	# seeds3 = catFq()
	# seeds = seeds3
	# print 'After cat'
	# for key in seeds.keys():
		# print key, seeds[key]
	mergeRunFq()
	dedup()
	trimFq()
	fq2fa()
	bowtie()
	samsort()
	sam2wig()
	mpileup()
	finalpileup()
	mergePairFq()
	dna2prot()
	bowtie_human()
	sam2wig_human()
	#bowtie_multi()
	#sam2wig_multi()
	mpileup_precise()
	wig2svg()
	
	f=open('pipeline_run.sh', 'w')
	#f.write('source mergePairFq.sh\nwait\n')
	f.write('source mergeRunFq.sh\nwait\n')
	f.write('source trim2ends.sh\nwait\n')
	f.write('source fq2fa.sh\nwait\n')
	f.write('source dedup.sh\nwait\n')
	f.write('source dna2prot.sh\nwait\n')
	f.write('source bowtie.sh\nwait\n')
	f.write('source sam2wig.sh > coverage.txt\nwait\n')
	f.write('source wig2svg.sh \nwait\n')
	# f.write('source bowtie_human.sh\nwait\n')
	# f.write('source sam2wig_human.sh > coverage_human.txt\nwait\n')
	#f.write('source bowtie_multi.sh\nwait\n')
	#f.write('source sam2wig_multi.sh > coverage_multi.txt\nwait\n')
	f.write('source samsort.sh\nwait\n')
	#f.write('source pileup.sh\nwait\n')
	f.write('source pileup_precise.sh\nwait\n')
	f.write('source finalpileup.sh\nwait\n')
	f.write('pileup2prot.py '+ref+' percent.pileup percent.prot.pileup\n')
	f.write('pileup2prot.py '+ref+' count.pileup count.prot.pileup\n')
	f.close()


