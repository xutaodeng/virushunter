#!/usr/bin/env python
from collections import defaultdict
import os.path
import os
import string
import sys
import itertools
import subprocess


servers=[  'bsidna4', \
 'bsidna5', 'bsidna6','bsidna7','bsidna8', \
  'bsidna9', 'bsidna10','bsidna11', \
  'bsidna14', 'bsidna15', 'bsidna19', 'bsidna20', \
  'bsidna23', 'bsidna24', 'bsidna25', 'bsidna26', \
'bsidna27','bsidna28','bsidna29', 'bsidna30'
]
nservers =len(servers)
def serverInfo():
	os.system('rm server.info')
	ps=[]
	for server in servers:
		p=subprocess.Popen('ssh '+server+' get_CPU.py '+server+' >> server.info\n', shell=True)
		ps.append(p)
	exit_codes = [p.wait() for p in ps]
	for code in exit_codes:
		if code !=0: print 'server error', code; sys.exit()
	print 'exit_codes', exit_codes
	SI={}
	f=open('server.info', 'r')
	for line in f:
		server, CPU, RAM = line.strip().split()
		SI[server]=(CPU, RAM)
	return SI
SI=serverInfo()
print SI

superservers=['bsidna35']

bowtiepath='/mnt/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2' 
blastxpath='/mnt/cluster/xdeng/tools/ncbi-blast-2.2.31+/bin/blastx'
blastnpath='/mnt/cluster/xdeng/tools/ncbi-blast-2.2.31+/bin/blastn'
dustmaskerpath='/mnt/cluster/xdeng/tools/ncbi-blast-2.2.31+/bin/dustmasker'
samtools='/mnt/cluster/xdeng/tools/samtools-1.2/samtools'


def readSeeds2(reAssemble='no'):
	global seeds
	stats=defaultdict()
	seedfile='fastq/samples.txt'
	f=open(seedfile, 'r')
	seeds=defaultdict(list)
	for line in f:
		key2=line.strip().rsplit('.', 1)[0]
		key1=line.strip().rsplit('_',4)[0]
		if len(key1)<len(key2):
			key=key1
		else:
			key=key2
		#key=line.strip().rsplit('_',2)[0]
		ss=line.strip()
		if line.strip().endswith('fasta'): ss = ss.replace('fasta', 'fastq')
		elif line.strip().endswith('sam'): ss = ss.replace('sam', 'fastq')
		elif line.strip().endswith('bam'): ss = ss.replace('bam', 'fastq')
		if len(seeds[key])<2: 
			seeds[key].append(ss)
			stats[key]=defaultdict(list)
		else:
			key=key+'.2'
			seeds[key].append(ss)
			stats[key]=defaultdict(list)
	# for line in f:
		# #key=line.strip().split('_')[0]
		# key1=line.strip().rsplit('.', 1)[0]
		# key2=line.strip().split('_')[0]
		# #print 'key1', key1, 'key2', key2
		# if len(key1)< len(key2):
			# key=key1
		# else:
			# key=key2
		# #print 'key', key
		# stats[key]=defaultdict(list)
	#print stats.keys()
	f.close()
	print seeds
	return stats, seeds


def bowtie(fafile):
	global pair
	f = open('bowtie.sh', 'w')
	f.write('/mnt/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build '+fafile+' '+wd+'/allbow\n')

	job=0
	bowindex= wd+'/allbow'
	for (key, fqfiles) in seeds.items():
		fastqfile1=wd+'/fastq/'+fqfiles[0]
		if pair: fastqfile2=wd+'/fastq/'+fqfiles[1]
		sam=wd+'/fastq/'+ key+'.sam'
		serverTag = 'ssh '+servers[job%nservers]
		# at least 30 bp perfect match -L 30, default is 20 
		if pair: f.write(serverTag+' '+bowtiepath+' --fast-local --quiet --local --no-hd --reorder -p 7 -L 30 -x '+bowindex+' -1 '+fastqfile1+' -2 '+fastqfile2+' -S '+sam+' &\n')
		else: f.write(serverTag+' '+bowtiepath+' --fast-local --quiet --local --no-hd --reorder -p 7 -L 30 -x '+bowindex+' -U '+fastqfile1+' -S '+sam+' &\n')
		job+=1
		if job%nservers==0: f.write('wait\n')
	f2 = open('sam2wig.sh', 'w')
	job=0
	countfiles=[]
	for (key, fqfiles) in seeds.items():
		sam=wd+'/fastq/'+ key+'.sam'
		wig=wd+'/fastq/'+ key+'.wig'
		serverTag = 'ssh '+servers[job%nservers]
		f2.write(serverTag+' '+dirscr+'sam2wig.py '+sam+' &\n')
		countfiles.append(wig)
		job+=1
		if job%nservers==0: f2.write('wait\n')
	f2.write('wait\n')
	f.close()
	f2.close()



if __name__ == "__main__":
	pair=False
	
	dirscr='/mnt/cluster/xdeng/script/'
	wd = os.path.abspath(os.path.dirname('.')).replace('san2', 'cluster2')
	serverInfo()
	readSeeds2()
	bowtie('fastq/spadeExt_contig.fasta')

