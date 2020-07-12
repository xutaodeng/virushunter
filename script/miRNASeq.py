#!/usr/bin/env python
from collections import defaultdict
import os.path
import os
import string
import sys
from os import walk
from os import listdir

servers=['bsidna3','bsidna4','bsidna5','bsidna6','bsidna7','bsidna8','bsidna9','bsidna10', 'bsidna11',\
'bsidna12','bsidna13','bsidna14','bsidna15','bsidna16','bsidna17','bsidna18','bsidna19','bsidna20',\
'bsidna21','bsidna22', 'bsidna23','bsidna24', 'bsidna25','bsidna26', 'bsidna27','bsidna28','bsidna29', 'bsidna30', \
'bsidna32', 'bsidna33', 'bsidna34', 'bsidna35', 'bsidna36']
genomefa= '/mnt/cluster/hg19/Sequence/Bowtie2Index/genome.fa'
mature= '/mnt/cluster/hg19/Annotation/SmallRNA/mature.fa'
hairpin = '/mnt/cluster/hg19/Annotation/SmallRNA/hairpin.fa'
bowtieIndex = ' /mnt/cluster/hg19/Sequence/BowtieIndex/genome '
#SOAPPATH="/mnt/cluster2/xdeng/tools/soap2.20release/"
# REFINDEX='/mnt/cluster2/xdeng/tools/Homo_sapiens/UCSC/hg19/Annotation/SmallRNA/combine-miRNA.fa.index'
# BLASTINDEX='/mnt/cluster2/xdeng/tools/Homo_sapiens/UCSC/hg19/Annotation/SmallRNA/combine-miRNA_blastdb'
#BLASTPATH='/mnt/cluster/tools/ncbi-blast-2.2.27+/bin/'

nservers=len(servers)

def readdir(seedfile='./samples.txt'):
	print 'hello'
	seeds={}
	# ls -1d */
	f=open(seedfile, 'r')
	for line in f:
		directory=line.strip()
		print directory
		files=[wd+'/fastq/'+directory+'/'+x for x in listdir('./fastq/'+directory) if (x.endswith('.fq') or x.endswith('.fastq'))]
	# for (dirpath, dirnames, filenames) in walk(mypath):
		# files=[wd+'/'+dirpath[2:]+'/'+x for x in filenames]
		seeds[directory]=files
	#print seeds
	f.close()
	print seeds.keys()
	print '--------------------------------------------------'
	for seed in seeds.keys():
		print seed, seeds[seed]
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
		print directory
		files=[wd+'/fastq/'+directory+'/'+x for x in listdir('./fastq/'+directory) if (x.endswith('.gz') or x.endswith('.fastq')) ]
		files1=[x for x in files if '_R1_' in x]
		files2=[x for x in files if '_R2_' in x]
		files1.sort()
		files2.sort()
		assert len(files1)==len(files2)
		for i in xrange(len(files1)):
			assert files1[i].replace('_R1_', '_R2_')==files2[i]
		seeds[directory]=(files1, files2)

	f.close()
	print seeds.keys()
	print '--------------------------------------------------'
	for seed in seeds.keys():
		print seed, seeds[seed]
	return seeds
	
# def prep_index():
	# f1=open(hairpin, 'r')
	# f2=open(mature, 'r')
	# of=open('/mnt/cluster2/xdeng/tools/Homo_sapiens/UCSC/hg19/Annotation/SmallRNA/combine-miRNA.fa', 'w')
	# id1=set()
	# id2=set()
	# for line in f1:
		# if line.strip().startswith('>'):
			# id1.add(line.strip().split()[0].lower())
			# of.write(line)
		# else:
			# of.write(line.replace('U', 'T'))

	# dup=False
	# dupcount=0
	# for line in f2:
		# if line.strip().startswith('>'):
			# id =line.strip().split()[0].lower()
			# if id in id1:
				# dup=True
				# dupcount+=1
			# else:
				# dup=False
		# if not dup:
			# if line.strip().startswith('>'):
				# of.write(line)
			# else:
				# of.write(line.replace('U', 'T'))
	# # print 'dupcount', dupcount

	# f1.close()
	# f2.close()
	# of.close()
	# cmd = SOAPPATH+'2bwt-builder /mnt/cluster2/xdeng/tools/Homo_sapiens/UCSC/hg19/Annotation/SmallRNA/combine-miRNA.fa'
	# cmd2= BLASTPATH+'makeblastdb -in /mnt/cluster2/xdeng/tools/Homo_sapiens/UCSC/hg19/Annotation/SmallRNA/combine-miRNA.fa -dbtype nucl -parse_seqids -out /mnt/cluster2/xdeng/tools/Homo_sapiens/UCSC/hg19/Annotation/SmallRNA/combine-miRNA_blastdb'
	# print cmd2
	# os.system(cmd2)

# def readSeeds(seedfile='samples.txt'):
	# #ls -1d */ > samples.txt
	# seeds=defaultdict(list)
	# f=open(seedfile, 'r')
	# for line in f:
		# key = line.strip()
		# fafile= key+'clean.fa'
		# seeds[key]= wd+'/'+fafile
	# f.close()
	# print seeds
	# return seeds

def uncompress():
	f = open('uncompress.sh', 'w') 
	job=0
	for (key, fafile) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		cmd = '"cd '+wd+'/fastq/ && gzip -d '+fafile[0]+'"'
		f.write(serverTag+' '+cmd+' &\n')
		job+=1
	f.close()

# def SOAP():
	# f = open('soap.sh', 'w') 
	# job=0
	# for (key, fafile) in seeds.items():
		# serverTag = 'ssh '+servers[job%nservers]
		# cmd = SOAPPATH+'soap -v 0 -r 2 -M 0 -a '+fafile+' -D '+REFINDEX+' -o '+fafile+'.soap'
		# f.write(serverTag+' '+cmd+' &\n')
		# job+=1
		# if job%nservers==0:  f.write('wait\n')
	# f.close()

# def blastn():
	# f = open('blastn.sh', 'w') 
	# job=0
	# for (key, fafile) in seeds.items():
		# serverTag = 'ssh '+servers[job%nservers]
		# cmd = BLASTPATH+'blastn  -query '+fafile+' -db '+BLASTINDEX+' -num_threads 8 -outfmt 6 -max_target_seqs 100000 -out '+fafile+'.blast'
		# f.write(serverTag+' '+cmd+' &\n')
		# job+=1
		# if job%nservers==0:  f.write('wait\n')
	# f.close()
def miRDeep2():
	f = open('miRDeep2.sh', 'w') 
	f2 = open('combineResults.sh', 'w')
	job=0
	csvs=[]
	for (key, fafile) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		serverTag=''
		fq=fafile[0].rstrip('.gz')
		#cmd0 = 'cd '+wd
		cmd1 ='	mapper.pl '+fq+' -e -h -j -k TGGAATTCTCGG  -l 18 -m -p '+bowtieIndex+' -s '+fq+'_collapsed.fa -t '+ \
			fq+'.arf -v'
		cmd2 =' quantifier.pl -p '+hairpin+' -m '+mature+' -r '+fq+'_collapsed.fa -t hsa -y '+key
		cmd = ' && '.join([ cmd1, cmd2])
		#f.write(serverTag+' '+cmd+' \n')
		f.write(serverTag+' '+cmd1+' && ' +cmd2+' \n')
		job+=1
		csvs.append('miRNAs_expressed_all_samples_'+key+'.csv')
	csvs.sort()
	f2.write('/mnt/cluster/xdeng/script/miRNAresult.py '+' '.join(csvs)+' '+'all.miRNA.read.count.csv all.miRNA.normalize.csv\n')
	f.close()
	f2.close()


def expression():
	files=[]
	for (key, fafile) in seeds.items():
		files.append(fafile+'.blast')
	e = {}
	allref=set()
	of = sys.stdout
	for file in files:
		f=open(file, 'r')
		if not e.has_key(file):
			e[file]=defaultdict(int)
		for line in f:
			parts = line.strip().split()
			ref= parts[1]
			e[file][ref]+=1
			allref.add(ref)
		f.close()
	allref = list(allref)
	print len(allref)
	of.write ('sample\t')
	for ref in allref:
		of.write(ref+'\t')
	of.write('\n')
	for file in files:
		of.write(file+'\t')
		for ref in allref:
			of.write(str(e[file][ref])+'\t')
		of.write('\n')

if __name__ == "__main__":
	# prep_index()
	wd = os.path.abspath(os.path.dirname('.'))
	if not wd.startswith('/mnt/'): wd='/mnt'+wd
	seeds=readdir()
	uncompress()
	miRDeep2()
	sf=open('pipeline_run.sh', 'w')
	sf.write('uncompress.sh\nwait\n') 
	sf.write('source miRDeep2.sh \nwait\n')
	sf.write('source combineResults.sh \n')
	sf.close()
	#expression()

