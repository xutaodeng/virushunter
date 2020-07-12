#!/usr/bin/env python
from collections import defaultdict
import os.path
import os
import string
import sys

# /mnt/cluster2/xdeng/tools/trinityrnaseq_r20140717/Trinity --seqType fq --JM 100G --left reads_1.fq  --right reads_2.fq --CPU 6

servers=['bsidna4','bsidna5','bsidna6','bsidna7','bsidna8','bsidna9','bsidna10',\
'bsidna11','bsidna12','bsidna13','bsidna14','bsidna15','bsidna17','bsidna18','bsidna19','bsidna20',\
'bsidna22','bsidna27','bsidna29', 'bsidna30']

servers=['bsidna7', 'bsidna8', 'bsidna5', 'bsidna6']
#servers=['bsidna11','bsidna12','bsidna9','bsidna10']
# servers=['bsidna32']

idbapath= '/mnt/cluster2/xdeng/tools/idba-1.1.1/bin/idba'
celerapath='/mnt/cluster2/xdeng/tools/wgs-8.1/Linux-amd64/bin/'
masurcapath='/mnt/cluster2/xdeng/tools/MaSuRCA-2.2.0/bin/'
omegapath='/mnt/cluster/tools/omega/omega/'
trinitypath='/mnt/cluster/xdeng/tools/trinityrnaseq-2.2.0/'
soappath='/mnt/cluster/tools/SOAPdenovo-63mer'
velvetg='/mnt/cluster/tools/velvet_1.2.10/velvetg'
velveth='/mnt/cluster/tools/velvet_1.2.10/velveth'
meta_velvetg="/mnt/cluster/tools/MetaVelvet-1.2.02/meta-velvetg"
Minimo="/mnt/cluster/tools/amos-3.1.0/bin/Minimo"
bowtiepath='/mnt/san/cluster/tools/bowtie2-2.1.0/bowtie2' 
blastxpath='/mnt/cluster/tools/ncbi-blast-2.2.27+/bin/blastx'
bowtieindexpath='/mnt/cluster/xdeng/hg19/mrnadna_bowtie'
abysspath='/mnt/cluster/tools/abyss/bin/abyss-pe'
bowtiebacs= ['/mnt/san/cluster/xdeng/nt/Bacteria1', \
			'/mnt/san/cluster/xdeng/nt/Bacteria2', \
			'/mnt/san/cluster/xdeng/nt/Bacteria3', \
			'/mnt/san/cluster/xdeng/nt/Bacteria4']
virusdbpath='/mnt/cluster/xdeng/blastdb/virus/virus_mask'
nvnrdbpath='/mnt/san/cluster/xdeng/blastdb/refseq/nvrefseq'
dirscr='/mnt/cluster/xdeng/script/'
Raytool='/mnt/cluster/tools/Ray2.3/Ray'
nservers=len(servers)


def genseedfile():
	os.system('cd fastq && ls -1 *fastq.gz > samples.txt && cd ..  ')

def readSeeds1(seedfile='fastq/samples.txt'):
	seeds=defaultdict(list)
	f=open(seedfile, 'r')
	for line in f:
		#key2=line.strip().rsplit('.', 1)[0]
		key2=line.strip().split('.')[0]
		key1=line.strip().split('_')[0]
		if len(key1)<len(key2):
			key=key1
		else:
			key=key2
		seeds[key].append(line.strip())
	f.close()
	return seeds

def readSeeds2(seedfile='fastq/samples.txt'):
	seeds=defaultdict(list)
	f=open(seedfile, 'r')
	for line in f:
		#key2=line.strip().rsplit('.', 1)[0]
		key2=''.join(line.strip().split('.')[0:2])
		key1=line.strip().split('_')[0]
		if len(key1)<len(key2):
			key=key1
		else:
			key=key2
		seeds[key].append(line.strip())
	f.close()
	return seeds

# def skipbowtie():
	# f = open('skipbowtie.sh', 'w')
	# for (key, fqfiles) in seeds.items():
		# i=1
		# for fqfile in fqfiles:
			# index=key+'_'+str(i)
			# fqfil=index+'.fil'
			# if fqfile.endswith('.gz'): f.write('zcat fastq/'+fqfile +' > fastq/'+fqfil+'\n')
			# else: f.write('cat fastq/'+fqfile +' > fastq/'+fqfil+'\n')
			# i+=1
	# f.close()
	
def skip_adaptor():
	f = open('skipadaptor.sh', 'w')
	for (key, fqfiles) in seeds.items():
		i=1
		for fqfile in fqfiles:
		#ff = open(wd+'/soap_config/'+key+'_soap.config', 'w')
			index=key+'_'+str(i)
			fid = wd+'/fastq/'+index
			f.write('recodeID.py '+fid+'.dup '+' '+fid+'.trim '+base+'_'+key+ ' '+str(i)+'\n') #quality trimming
			i+=1
	f.close()

def sampleFastq():
	f = open('sampleFastq.sh', 'w')
	for (key, fqfiles) in seeds.items():
		i=1
		fqs=[]
		for fqfile in fqfiles:
			fqs.append(fqfile)
			i+=1
		f.write('sampleFastq.py '+wd+'/fastq/'+fqs[0]+' '+wd+'/fastq/'+fqs[1]+' '+wd+'/fastq/out'+fqs[0][:-3]+' '+wd+'/fastq/out'+fqs[1][:-3]+' 75000\n') #quality trimming
	f.close()

def prepBlastFile_adaptor(): #for adaptor filtering
	try: os.mkdir(wd+'/blast_filter_out/')
	except: pass
	try: os.mkdir(wd+'/blast_filter_out/blast2/')
	except: pass
	f = open('prepBlastFile_adaptor.sh', 'w')
	path='fastq/'

	bpath=wd+'/blast_filter_out/blast2/'
	directory = os.path.basename(wd)
	alldb=[]
	job=0
	for (key, fqfiles) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		job+=1
		i=1
		pairdb=[]
		f.write(serverTag+' "cd '+wd+'/fastq/ &&')
		for fqfile in fqfiles:
			index=key+'_'+str(i)
			fqfil=index+'.fil'
			f.write('makeblastdb -dbtype nucl -parse_seqids -in '+index+'.fa -out '+index+'_blastdb && ')
			pairdb.append(index+'_blastdb')
			i+=1
		f.write('blastdb_aliastool -dblist '+'\''+' '.join(pairdb)+'\' -dbtype nucl -out '+key+'_blastdb -title \''+key+'_blastdb\' "&\n')
		if job%nservers==0:  f.write('wait\n')
		alldb.append(key+'_blastdb')
	f.write('wait\n')
	f.write('blastdb_aliastool -dblist '+'\''+' '.join(alldb)+'\' -dbtype nucl -out '+directory+'_blastdb -title \''+directory+'_blastdb\' \n')
	f.write('mv fastq/*_blastdb* '+bpath+'\n')
	#f.write('cd ..\n')
	f.close()

def trim():
	f = open('clonetrim.sh', 'w')
	f1 = open('clonefq2fa.sh', 'w')
	f2 = open('blast_adaptor.sh', 'w')
	f3 = open('blasttrim.sh', 'w')
	f4 = open('qualitytrim.sh', 'w')
	f5 = open('fq_clean.sh', 'w')
	adaptors = dirscr+'adaptor.fa'
	bpath=wd+'/blast_filter_out/blast2/'
	job=0
	for (key, fqfiles) in seeds.items():
		i=1
		trimfiles=[]
		seqfiles=[]
		for fqfile in fqfiles:
		#ff = open(wd+'/soap_config/'+key+'_soap.config', 'w')
			index=key+'_'+str(i)
			fid = wd+'/fastq/'+index
			fastqfile=wd+'/fastq/'+fqfile
			fqfil=fid+'.fil'
			serverTag = 'ssh '+servers[job%nservers]
			job+=1
			f.write(serverTag + ' '+dirscr+'dedup.py '+fqfil+' '+fid+'.dup &\n')
			f1.write(dirscr+'fq2faID.py '+fid+'.dup '+' '+index+' '+fid+'.fa\n')
			f2.write(serverTag + ' blastn -task blastn -evalue 1  -max_target_seqs 100000000 -outfmt \'"6  qseqid  sseqid evalue qstart qend sstart send"\'')
			f2.write(' -query '+adaptors+ ' -num_threads '+thread+'  -db '+bpath+index+'_blastdb'+' -out '+fid+'.tab & \n')
			f3.write(serverTag + ' '+dirscr+'blast_trim.py '+fid+'.dup '+fid+'.tab '+fid+'.ada '+ base+'_'+key+ ' '+str(i)+' &\n') #adaptor trimming
			#f4.write(serverTag + ' '+dirscr+'trim_quality.py '+fid+'.ada '+fid+'_sequence.txt 33 '+fastqfile+' &\n') #quality trimming
			f4.write(serverTag + ' '+dirscr+'trim_quality.py '+fid+'.ada '+fid+'.trim 33 '+fastqfile+' &\n') #quality trimming
			trimfiles.append(fid+'.trim')
			seqfiles.append(fid+'_sequence.txt')
			if job%nservers==0: 
				f.write('wait\n')
				f1.write('wait\n')
				f2.write('wait\n')
				f3.write('wait\n')
				f4.write('wait\n')
				f5.write('wait\n')
			i+=1
		f5.write(serverTag + ' '+dirscr+'fq_pair_clean.py '+' '.join(trimfiles)+' '+' '.join(seqfiles)+' 5 &\n')
	f.close()
	f1.close()
	f2.close()
	f3.close()
	f4.close()
	f5.close()

def check_pair_overlap():
	f = open('check_pair.sh', 'w')
	job=0
	for (key, fqfiles) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		job+=1
		f.write(serverTag + ' '+dirscr+'check_pair.py '+wd+'/fastq/' + key+'_1_sequence.txt '+wd+'/fastq/' + key+'_2_sequence.txt '+' &\n') #quality trimming
		if job%nservers==0: 
			f.write('wait\n')
	f.close()

def fq_check():
	f = open('fq_check.sh', 'w')
	job=0
	for (key, fqfiles) in seeds.items():
		for fqfile in fqfiles:
			serverTag = 'ssh '+servers[job%nservers]
			job+=1
			f.write(serverTag + ' '+dirscr+'fq_check.py '+wd+'/fastq/' + fqfile+' &\n') #quality trimming
			if job%nservers==0: 
				f.write('wait\n')
	f.close()
	
def clean_dir():
	f2 = open(wd+'/clean.sh', 'w')
	#f2.write('find ./fastq/ -type f ! -name "*.gz"  -a ! -name "*.fastq" ! -name "*.log" ! -name "*.php" ! -name "*.gif" -delete\nwait\n')
	f2.write('rm -rf blast_* \nwait\n')
	f2.write('rm -rf soap_* \nwait\n')
	f2.write('rm -rf *.time \nwait\n')
	f2.write('rm -rf velvet_* \nwait\n')
	f2.write('rm -rf run \nwait\n')
	f2.write('rm -rf *.sh \nwait\n')
	f2.write('rm -rf *.log \nwait\n')
	f2.close()

def polyA():
	job=0
	f = open('polyA_raw.sh', 'w')
	f2 = open('polyA_clean.sh', 'w')
	try: os.mkdir(wd+'/'+base)
	except: pass
	try: os.mkdir(wd+'/'+base+'/hist/')
	except: pass
	for (key, fqfiles) in seeds.items():
		i=1
		for fqfile in fqfiles:
			serverTag = 'ssh '+servers[job%nservers]
			keypath=wd+'/'+base+'/hist/'+key
			f.write(serverTag + ' '+ dirscr+'polyA.py '+wd+'/fastq/'+fqfile+' '+keypath+'_'+str(i)+' '+'raw  &\n')
			f2.write(serverTag + ' '+dirscr+'polyA.py '+wd+'/fastq/'+ key+'_'+str(i)+'_sequence.txt '+keypath+'_'+str(i)+' clean &\n')
			job+=1
			if job%nservers==0: 
				f.write('wait\n')
				f2.write('wait\n')
			i+=1
	f.close()
	f2.close()

def skipbowtie():
	f = open('skipbowtie.sh', 'w')
	for (key, fqfiles) in seeds.items():
		i=1
		for fqfile in fqfiles:
			index=key+'_'+str(i)
			fqfil=index+'.fil'
			if fqfile.endswith('.gz'): f.write('zcat fastq/'+fqfile +' > fastq/'+fqfil+'\n')
			else: f.write('cat fastq/'+fqfile +' > fastq/'+fqfil+'\n')
			i+=1
	f.close()

# def bowtieNT():
	# f = open('bowtieNT.sh', 'w')
	# f2 = open('sam2fq.sh', 'w')
	# #f2.write('rm stat.log\n')
	# job=0
	# for (key, fqfiles) in seeds.items():
		# i=1
		# for fqfile in fqfiles:
			# index=key+'_'+str(i)
			# for j in range(1, 13):
				# serverTag = 'ssh '+servers[job%nservers]
				# f.write(serverTag+' '+bowtiepath+' --no-hd --quiet --local -p 8 -x /mnt/cluster/xdeng/nt/nt'+str(j)+ \
				# ' -U '+wd+'/fastq/'+index +'_sequence.txt.tmp -S '+wd+'/fastq/'+ index+'_'+str(j)+'.sam &\n')
				# job+=1
				# if job%nservers==0: f.write('wait\n')
			# f2.write(dirscr+'sam2fq.py '+wd+'/fastq/' + index+' '+wd+'/fastq/'+index+'_sequence.txt \n')
			# i+=1
	# f.close()
	# f2.close()

def sam2fq(pair, keep_human):
	f2 = open('bowtiesam2fq.sh', 'w')
	job2=0
	for (key, fqfiles) in seeds.items():
		fqfils=[]
		i=1
		for fqfile in fqfiles:
			index=key+'_'+str(i)
			fqfil=wd+'/fastq/'+index+'.fil'
			fqfils.append(fqfil)
			i+=1
		serverTag2 = 'ssh '+servers[job2%nservers]
		if not keep_human : start, end = '0', '4' #remove human
		if keep_human : start, end = '1', '4'
		if pair: f2.write(serverTag2 + ' '+dirscr+'sam2fq_bac.py '+wd+'/fastq/'+key+' '+start+' '+end+' '+fqfils[0]+' '+fqfils[1]+ ' & \n')
		else: f2.write(serverTag2 + ' '+dirscr+'sam2fq_bac.py '+wd+'/fastq/'+key+' '+start+' '+end+' '+fqfils[0]+ ' & \n')
		job2+=1
		if job2%nservers==0: f2.write('wait\n')
	f2.close()

def bowtieBac(pair, keep_human):
	f = open('bowtieBac.sh', 'w')
	job=0
	if keep_human: bowindex=bowtiebacs
	else: bowindex= [bowtieindexpath]; bowindex.extend(bowtiebacs)
	print bowindex
	for (key, fqfiles) in seeds.items():
		if keep_human: bj=1 #start with bacteria 1,2,3,4
		else: bj=0 #remove human start with 0
		for bowtieind in bowindex:
			i=1
			for fqfile in fqfiles:
				fastqfile=wd+'/fastq/'+fqfile
				index=key+'_'+str(i)
				sam=wd+'/fastq/'+ index+'_'+str(bj)+'.sam'
				serverTag = 'ssh '+servers[job%nservers]
				f.write(serverTag+' '+bowtiepath+' --quiet --local --no-hd --reorder -p 7 -x '+bowtieind+' -U '+fastqfile+' -S '+sam+' &\n')
				job+=1
				if job%nservers==0: f.write('wait\n')
				i+=1
			bj+=1
	#f.write('exit')
	f.close()

def bowtieRef(blastRef):
	f = open('bowtieRef.sh', 'w')
	if blastRef== 'NIBSC': 
		bowindex= '/mnt/cluster/xdeng/blastdb/NIBSC'
	elif blastRef== 'ref5': 
		bowindex= '/mnt/cluster/xdeng/blastdb/ref5Virus'
	elif blastRef== 'BASV': 
		bowindex= '/mnt/cluster/xdeng/blastdb/BASV'
	elif blastRef== 'bacteria': 
		bowindex= '/mnt/cluster/xdeng/blastdb/bacteria'
	elif blastRef== 'SURPI': 
		bowindex= '/mnt/cluster/xdeng/blastdb/SURPI'
	for (key, fqfiles) in seeds.items():
		i=1
		sams=[]
		for fqfile in fqfiles:
			fastqfile=wd+'/fastq/'+fqfile
			fastqfile=wd+'/fastq/'+key+'_'+str(i)+'_sequence.txt'
			index=key+'_'+str(i)
			sam=wd+'/fastq/'+ index+'_ref.sam'
			f.write(bowtiepath+' --quiet --local --no-hd --reorder -p 7 -x '+bowindex+' -U '+fastqfile+' -S '+sam+' \n')
			sams.append(sam)
			i+=1
		samfile=wd+'/fastq/'+ key+'_ref.sam'
		wigfile=wd+'/fastq/'+ key+'_ref.wig'
		svgfile=wd+'/fastq/'+ key+'_ref.svg.html'  
		f.write('cat '+' '.join(sams)+ ' > '+samfile+'\n')
		# f.write(dirscr+'sam2fq_BASV.py '+sams[0]+' '+sams[1]+' \n')
		f.write(dirscr+'sam2wig.py '+samfile+' '+wigfile+' \n')
		f.write(dirscr+'wig2svg.py '+wigfile+' '+svgfile+' 50 \n') #resolution
	f.close()
	
def prep_sra():
	f = open('sra.sh', 'w')
	for (key, fqfiles) in seeds.items():
		i=1
		sams=[]
		fqfils=[]
		for fqfile in fqfiles:
			fastqfile=wd+'/fastq/'+fqfile
			f.write(dirscr+'sra.py '+fastqfile+' '+key+'_'+str(i)+' '+wd+'/sra.fq.gz\n')
			i+=1
	f.close()



def prep_reads(n, length, pair): #filter length fq file and contig, and then combine
	f = open(wd+'/prep_reads.sh', 'w')
	job=0
	if skipall:
		for (key, fqfiles) in seeds.items(): # copy sequence.txt to fq for Ray input requirement
			serverTag = 'ssh '+servers[job%nservers]
			f.write(serverTag+' cat '+wd+'/fastq/'+fqfiles[0]+ ' > '+wd+'/fastq/'+key+'_1_sequence.txt &\n')
			if pair: f.write(serverTag+' cat '+wd+'/fastq/'+fqfiles[1]+ ' > '+wd+'/fastq/'+key+'_2_sequence.txt &\n')
			job+=1
			if job%nservers==0: f.write('wait\n')
		f.write('wait\n')
	#merge .fq
	job=0
	for (key, fqfiles) in seeds.items():
		file1 = wd+'/fastq/'+key+'_1_sequence.txt'
		if pair: file2 =wd+'/fastq/'+key+'_2_sequence.txt'
		serverTag = 'ssh '+servers[job%nservers]
		f.write(serverTag+' cat '+ file1+' ')
		if pair: f.write(file2+' ')
		f.write('> '+wd+'/fastq/'+key+'.fq &\n')
		job+=1
		if job%nservers==0: f.write('wait\n')
	f.write('wait\n')

	#prep .fa
	job=0
	for (key, vdict) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		f.write(serverTag+' '+dirscr+'fq2fa.py ')
		f.write(wd+'/fastq/'+key+'.fq '+wd+'/fastq/'+key+'.fa '+str(length)+' &\n')
		job+=1
		if job%nservers==0: f.write('wait\n')
	f.write('wait\n')
	#below prep abyss
	job=0
	for (key, vdict) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		f.write(serverTag+' '+dirscr+'fqLenFilter.py '+ wd+'/fastq/'+key+'.fq '+' '+wd+'/fastq/'+key+'abyss.fq 35\n')
		job+=1
		if job%nservers==0: f.write('wait\n')
	#below prep raw.fq, untrimed raw fq
	# job=0
	# for (key, fqfiles) in seeds.items():
		# serverTag = 'ssh '+servers[job%nservers]
		# i=1
		# for fqfile in fqfiles:
			# index=key+'_'+str(i)
			# raw=wd+'/fastq/'+index+'raw.fq'
			# if fqfile.endswith('.gz'): f.write(serverTag+' zcat '+wd+'/fastq/'+fqfile +' > '+raw+' &\n')
			# else: f.write(serverTag+' cat '+wd+'/fastq/'+fqfile +' > '+raw+' &\n')
			# i+=1
			# job+=1
			# if job%nservers==0: f.write('wait\n')
	f.write('wait\n')
	#below prep mira, change labels because they are too long for mira to handle
	job=0
	for (key, fqfiles) in seeds.items():
		i=1
		for fqfile in fqfiles:
			serverTag = 'ssh '+servers[job%nservers]
			index=key+'_'+str(i)
			# raw=wd+'/fastq/'+index+'raw.fq'
			seqFile = wd+'/fastq/'+index+'_sequence.txt'
			fqout=wd+'/fastq/'+index+'mira.fastq'
			f.write(serverTag+' '+dirscr+'mira_shortenID.py '+ seqFile+' '+fqout+' '+str(i)+':N:0:12 &\n')
			#f.write(serverTag+' '+dirscr+'mira_shortenID.py '+ raw+' '+fqout+' &\n')
			i+=1
			job+=1
			if job%nservers==0: f.write('wait\n')
	f.close()


def meta_velvet():
	f = open('meta_velvet.sh', 'w')
	# velveth out-dir 51 -fastq -shortPaired HMP.small/SRR041654_shuffled.fastq  HMP.small/SRR041655_shuffled.fastq 
	# velvetg out-dir -exp_cov auto -ins_length 260
	# meta-velvetg out-dir -ins_length 260 | tee logfile
	rval=[]
	job=0
	for (key, vdict) in seeds.items():
		outdir = wd+'/velvet_'+key+'/'
		serverTag = '{ time ssh '+servers[job%nservers]
		#serverTag=''
		if pair: 
			f.write(serverTag+' "'+velveth+' '+outdir+' '+metakmer)
			f.write(' -shortPaired -fastq '+wd+'/fastq/'+key+'_1_sequence.txt '+wd+'/fastq/'+key+'_2_sequence.txt && ')
			f.write(velvetg+' '+outdir+' -exp_cov auto -ins_length 260 && ')
			f.write(meta_velvetg+' '+outdir+' -ins_length 260"  ; } 2> '+wd+'/'+key+'_meta.time &\n')
		else:
			f.write(serverTag+' "'+velveth+' '+outdir+' '+metakmer)
			f.write(' -short -fastq '+wd+'/fastq/'+key+'_1_sequence.txt && ')
			f.write(velvetg+' '+outdir+' -exp_cov auto && ')
			f.write(meta_velvetg+' '+outdir+'"  ; } 2> '+wd+'/'+key+'_velvet.time &\n')
		velvetOut=wd+'/velvet_'+key+'/contigs.fa'
		rval.append(velvetOut)
		job+=1
		if job%nservers==0: f.write('wait\n')
	f.close()
	return rval

def trinity(RAM):
	f = open('trinity.sh', 'w')
	rval=[]
	job=0
	for (key, vdict) in seeds.items():
		outdir = wd+'/trinity_'+key+'/'
		serverTag = '{ time ssh '+servers[job%nservers]
		file1 = wd+'/fastq/'+key+'_1_sequence.txt'
		file2 = wd+'/fastq/'+key+'_2_sequence.txt'
		
		if pair: 
			#Trinity --seqType fq --JM 100G --left reads_1.fq  --right reads_2.fq --CPU 6
			f.write(serverTag+' "'+trinitypath+'/Trinity --seqType fq --JM '+RAM+' --output '+outdir+' --left '+file1 + \
					' --right '+file2+' --CPU 8 "  ; } 2> '+wd+'/'+key+'_trinity.time &\n')
		else:
			f.write(serverTag+' "'+trinitypath+'/Trinity --seqType fq --JM '+RAM+' --output '+outdir+' --single '+file1 + \
					' --CPU 8 "  ; } 2> '+wd+'/'+key+'_trinity.time &\n')
		trinityOut=wd+'/trinity_'+key+'/Trinity.fasta'
		rval.append(trinityOut)
		job+=1
		if job%nservers==0: f.write('wait\n')
	f.close()
	return rval
	
def masurca():
	f = open('masurca.sh', 'w')
	rval=[]
	job=0
	for (key, vdict) in seeds.items():
		outdir = wd+'/masurca_'+key+'/'
		try: os.mkdir(outdir)
		except: pass
		serverTag = '{ time ssh '+servers[job%nservers]
		file1 = wd+'/fastq/'+key+'_1_sequence.txt'
		file2 = wd+'/fastq/'+key+'_2_sequence.txt'
		f1 = open(outdir+'config.txt', 'w')
		f1.write('DATA\n')
		f1.write('PE= pe 180 20  '+file1+' '+file2+'\n')
		f1.write('END\n')
		f1.write('PARAMETERS\n')
		f1.write('GRAPH_KMER_SIZE=auto\n')
		f1.write('USE_LINKING_MATES=1\n')
		f1.write('LIMIT_JUMP_COVERAGE = 60\n')
		f1.write('CA_PARAMETERS = ovlMerSize=30 cgwErrorRate=0.25 ovlMemory=4GB\n')
		f1.write('KMER_COUNT_THRESHOLD = 1\n')
		f1.write('NUM_THREADS= 8\n')
		f1.write('JF_SIZE=100000000\n')
		f1.write('DO_HOMOPOLYMER_TRIM=0\n')
		f1.write('END\n')
		f1.close()
		f.write(serverTag+' "cd '+ outdir + ' && ' + masurcapath+'masurca config.txt && '+outdir+'assemble.sh"  ; } 2> '+wd+'/'+key+'_masurca.time &\n')
		masurcaOut=outdir+'CA/10-gapclose/genome.ctg.fasta'
		rval.append(masurcaOut)
		job+=1
		if job%nservers==0: f.write('wait\n')
	f.close()
	return rval


	
def celera(): #filter length fq file and contig, and then combine
	#/mnt/cluster2/xdeng/tools/wgs-8.1/Linux-amd64/bin> ./fastqToCA -libraryname TEST -technology none -reads test.fastq  > test.frg
	#./runCA -d test -p test test.frg
	f0 = open(wd+'/celera_prepReads.sh', 'w')
	f = open(wd+'/celeraPrep.sh', 'w')
	f1 = open(wd+'/celera.sh', 'w')
	f2 = open(wd+'/celeraCombine.sh', 'w')
	rval=[]
	job=0
	for (key, fqfiles) in seeds.items():
		i=1
		#rawfq = wd+'/fastq/'+key+'.raw.fq'
		rawfq = wd+'/fastq/'+key+'abyss.fq'
		fqfiles2 = [wd+'/fastq/'+fqfile for fqfile in fqfiles]
		serverTag = '{ time ssh '+servers[job%nservers]
		serverTag0 = '{ time '
		f0.write('cat '+' '.join(fqfiles2)+' > '+rawfq+' & \n')
		label=key+'CELERA'
		f.write('cd '+wd+' && '+celerapath+'fastqToCA  -libraryname '+ label +' -technology none -reads '+rawfq+' > '+label+'.frg\n')
		#f1.write(serverTag+' "cd '+wd+' && '+celerapath+'runCA -d '+label+' -p '+label+' '+label+'.frg" ; } 2> '+wd+'/'+key+'_celera.time &\n')
		f1.write(serverTag0+' '+celerapath+'runCA -d '+label+' -p '+label+' '+label+'.frg ; } 2> '+wd+'/'+key+'_celera.time \n')
		prefix = label+'/9-terminator/'+label
		fasta = [prefix+'.ctg.fasta', prefix+'.deg.fasta', prefix+'.utg.fasta', prefix+'.scf.fasta']
		combinefasta = label+'/celera.contig.fa'
		f2.write('cd '+wd+' && cat '+' '.join(fasta)+' > '+combinefasta+'\n')
		rval.append(wd+'/'+combinefasta)
		job+=1
		if job%nservers==0: f1.write('wait\n')
	f.close()
	f1.close()
	f2.close()
	f0.close()
	return rval
# now generate commandline sh for SOAPDENOVO
def soap_single():
	f = open('soap.sh', 'w')
	try: os.mkdir(wd+'/soap_config/')
	except: pass
	try: os.mkdir(wd+'/soap_out/')
	except: pass
	job=0
	rval=[]
	for (key, vdict) in seeds.items():
		serverTag = '{ time ssh '+servers[job%nservers]
		ff = open(wd+'/soap_config/'+key+'_soap.config', 'w')
		f.write(serverTag+' "'+soappath+' all -K '+ soapkmer)
		f.write(' -s '+ wd+'/soap_config/'+key+'_soap.config ')
		f.write(' -R -o '+ wd+'/soap_out/'+key+'_soap"  ; } 2> '+wd+'/'+key+'_soap.time &\n')
		ff.write('max_rd_len=300\n[LIB]\n')
		ff.write('nreverse_seq=0\nasm_flags=3\nrank=1\n')
		ff.write('q='+wd+'/fastq/'+key+'_1_sequence.txt\n')
		ff.close()
		soapOut=wd+'/soap_out/'+key+'_soap.contig'
		rval.append(soapOut)
		job+=1
		if job%nservers==0: f.write('wait\n')
	f.close()
	return rval
	
def soap_pair():
	f = open('soap.sh', 'w')
	try: os.mkdir(wd+'/soap_config/')
	except: pass
	try: os.mkdir(wd+'/soap_out/')
	except: pass
	job=0
	rval=[]
	for (key, vdict) in seeds.items():
		serverTag = '{ time ssh '+servers[job%nservers]
		ff = open(wd+'/soap_config/'+key+'_soap.config', 'w')
		f.write(serverTag+' "'+soappath+' all -K '+soapkmer)
		f.write(' -s '+ wd+'/soap_config/'+key+'_soap.config ')
		f.write(' -R -o '+ wd+'/soap_out/'+key+'_soap"  ; } 2> '+wd+'/'+key+'_soap.time &\n')
		ff.write('max_rd_len=300\n[LIB]\n')
		ff.write('avg_ins=350\nreverse_seq=0\nasm_flags=3\nrank=1\n')
		ff.write('q1='+wd+'/fastq/'+key+'_1_sequence.txt\n'+'q2='+wd+'/fastq/'+key+'_2_sequence.txt \n')
		ff.close()
		soapOut=wd+'/soap_out/'+key+'_soap.contig'
		rval.append(soapOut)
		job+=1
		if job%nservers==0: f.write('wait\n')
	f.close()
	return rval

def omega(): #filter length fq file and contig, and then combine
	f = open(wd+'/omega.sh', 'w')
	cat=[]
	rval=[]
	job=0
	for (key, vdict) in seeds.items():
		try: os.mkdir(wd+'/abyss_'+key)
		except: pass
		serverTag = '{ time ssh '+servers[job%nservers]
		serverTag0 = '{ time '
		# /mnt/cluster/tools/omega/omega -se 1 /mnt/cluster2/xdeng/Vhunt/140527_ZW_plasma/fastq/pool4abyss.fq -f test -l 50
		f.write(serverTag+' "cd '+wd+' && /mnt/cluster/tools/omega/omega -se 1 '+wd+'/fastq/'+key+'abyss.fq -l 60 -f '+key+'omega" ; } 2> '+wd+'/'+key+'_omega.time &\n')
		# f.write(serverTag+' "cd '+wd+' && abyss-pe -C '+wd+'/abyss_'+key+' name='+key+' k='+abysskmer+' se='+wd+'/fastq/'+key+'abyss.fq" ; } 2> '+wd+'/'+key+'_abyss.time &\n')
		cat.append('cat '+key+'omegacontigs*.fasta > '+key+'-omega.fa' )
		rval.append(wd+'/'+key+'-omega.fa')
		job+=1
		if job%nservers==0: f.write('wait\n')
	f.write('wait\n'+'\n'.join(cat)+'\n')
	f.close()
	return rval

def abyss(): #filter length fq file and contig, and then combine
	f = open(wd+'/abyss.sh', 'w')
	rval=[]
	job=0
	for (key, vdict) in seeds.items():
		try: os.mkdir(wd+'/abyss_'+key)
		except: pass
		serverTag = '{ time ssh '+servers[job%nservers]
		f.write(serverTag+' "cd '+wd+' && abyss-pe -C '+wd+'/abyss_'+key+' name='+key+' k='+abysskmer+' se='+wd+'/fastq/'+key+'abyss.fq" ; } 2> '+wd+'/'+key+'_abyss.time &\n')
		rval.append(wd+'/abyss_'+key+'/'+key+'-unitigs.fa')
		job+=1
		if job%nservers==0: f.write('wait\n')
	f.close()
	return rval

def idba(): #filter length fq file and contig, and then combine
	##/mnt/cluster2/xdeng/tools/idba-1.1.1/bin/idba -r SD9.fa -o hello.fa
	f = open(wd+'/idba.sh', 'w')
	rval=[]
	job=0
	for (key, vdict) in seeds.items():
		serverTag = '{ time ssh '+servers[job%nservers]
		cmd = idbapath+' -r '+wd+'/fastq/'+key+'.fa -o '+wd+'/'+key+'_idba'
		f.write(serverTag+' "'+cmd+'" ; } 2> '+wd+'/'+key+'_idba.time &\n')
		rval.append(wd+'/'+key+'_idba/contig.fa')
		job+=1
		if job%nservers==0: f.write('wait\n')
	f.close()
	return rval
# def Ray(): #filter length fq file and contig, and then combine
	# f = open(wd+'/Ray.sh', 'w')
	# rval=[]
	# for (key, vdict) in seeds.items():
		# fq1 = wd+'/fastq/'+key+'_1.fq '
		# fq2 = wd+'/fastq/'+key+'_2.fq '
		# f.write('mpiexec -n 1 '+Raytool+' -o '+key+'_ray -k 31 -p '+fq1+' '+fq2+'\n')
		# rval.append(wd+'/'+key+'-unitigs.fa')
	# f.close()
	# return rval

def partition(): #filter length fq file and contig, and then combine
	f = open(wd+'/partition.sh', 'w')
	for (key, vdict) in seeds.items():
		f.write(dirscr+'partition.py '+wd+'/fastq/'+key+'abyss.fq  100000\n')
	f.close()

def abyss_partition(): #filter length fq file and contig, and then combine
	f = open(wd+'/abyss_partition.sh', 'w')
	f1 = open(wd+'/abyss_combine.sh', 'w')
	f2=open(wd+'/abyss_partition_blast.sh', 'w')
	rval=[]
	job =0
	for (key, vdict) in seeds.items():
		rv=[]
		for j in xrange(20): #assume 100 chunks, unknown file size but some may not exist
			try: os.mkdir(wd+'/abyss_'+key+'_chunk'+str(j))
			except: pass
			serverTag = 'ssh '+servers[job%nservers]
			f.write(serverTag+' "cd '+wd+' && abyss-pe -C '+wd+'/abyss_'+key+'_chunk'+str(j)+' name='+key+'_chunk'+str(j)+ \
			' k='+abysskmer+' se='+wd+'/fastq/'+key+'abyss.fq_'+str(j)+'" &\n')
			#f.write('abyss-pe name='+key+'_chunk'+str(j)+' k='+abysskmer+' se='+wd+'/fastq/'+key+'.fq_'+str(j)+'  \n')
			outfile = wd+'/abyss_'+key+'_chunk'+str(j)+'/'+key+'_chunk'+str(j)+'-unitigs.fa'
			rv.append(outfile)
			if 'setC' in key and j<=5:
				f2.write(dirscr+'blastContig.py '+outfile+'  abyss_chunk_'+str(j)+' '+blastRef+'\n')
			#rval.append(wd+'/abyss_'+key+'/'+key+'-unitigs.fa')
			job+=1
			if job%nservers==0: f.write('wait\n')
		f1.write('cat '+' '.join(rv) + ' > '+wd+'/'+key+'_abyss_partition.fa\n')
		rval.append(wd+'/'+key+'_abyss_partition.fa')
	f.close()
	f1.close()
	f2.close()
	return rval

def soap_partition(): #single end soap partition
	try: os.mkdir(wd+'/soap_partition_config/')
	except: pass
	try: os.mkdir(wd+'/soap_partition_out/')
	except: pass
	job=0
	rval=[]
	f = open(wd+'/soap_partition.sh', 'w')
	f1 = open(wd+'/soap_combine.sh', 'w')
	rval=[]
	for (key, vdict) in seeds.items():
		rv=[]
		for j in xrange(20): #assume 100 chunks, unknown file size but some may not exist
			serverTag = 'ssh '+servers[job%nservers]
			f.write(serverTag+' "'+soappath+' all -K '+ soapkmer)
			f.write(' -s '+ wd+'/soap_partition_config/'+key+'_'+str(j)+'_soap.config ')
			f.write(' -R -o '+ wd+'/soap_partition_out/'+key+'_'+str(j)+'_soap" &\n')
			ff = open(wd+'/soap_partition_config/'+key+'_'+str(j)+'_soap.config', 'w')
			ff.write('max_rd_len=300\n[LIB]\n')
			ff.write('nreverse_seq=0\nasm_flags=3\nrank=1\n')
			ff.write('q='+wd+'/fastq/'+key+'abyss.fq_'+str(j)+'  \n')
			ff.close()
			soapOut=wd+'/soap_partition_out/'+key+'_'+str(j)+'_soap.contig'
			rv.append(soapOut)
			job+=1
			if job%nservers==0: f.write('wait\n')
		f1.write('cat '+' '.join(rv) + ' > '+wd+'/'+key+'_soap_partition.fa\n')
		rval.append(wd+'/'+key+'_soap_partition.fa')
	f.close()
	f1.close()
	return rval

def velvet_partition(): #single end soap partition
	job=0
	rval=[]
	f = open(wd+'/velvet_partition.sh', 'w')
	f1 = open(wd+'/velvet_combine.sh', 'w')
	rval=[]
	for (key, vdict) in seeds.items():
		rv=[]
		for j in xrange(20): #assume 100 chunks, unknown file size but some may not exist
			outdir = wd+'/velvet_'+key+'_'+str(j)+'/'
			serverTag = 'ssh '+servers[job%nservers]
			f.write(serverTag+' "'+velveth+' '+outdir+' '+metakmer)
			f.write(' -short -fastq '+wd+'/fastq/'+key+'abyss.fq_'+str(j)+' && ')
			f.write(velvetg+' '+outdir+' -exp_cov auto && ')
			f.write(meta_velvetg+' '+outdir+'" &\n')
			velvetOut=wd+'/velvet_'+key+'_'+str(j)+'/contigs.fa'
			rv.append(velvetOut)
			job+=1
			if job%nservers==0: f.write('wait\n')
		f1.write('cat '+' '.join(rv) + ' > '+wd+'/'+key+'_velvet_partition.fa\n')
		rval.append(wd+'/'+key+'_velvet_partition.fa')
	f.close()
	f1.close()
	return rval

def mira4(mira_local):
	f = open('mira.sh', 'w')
	job=0
	rval=[]
	for (key, vdict) in seeds.items():
		serverTag = '{ time ssh '+servers[job%nservers]
		ff = open(wd+'/'+key+'.conf', 'w')
		f.write(serverTag+' "cd '+wd+' && mira -t 8 '+wd+'/'+key+'.conf" ; } 2> '+wd+'/'+key+'_mira.time &\n')
		ff.write('project = '+key+'\n')
		if mira_local: ff.write('parameters = -GE:not=8 -DI:trt=/home/BSI/306307/ -OUTPUT:rtd=yes\n')
		else: ff.write('parameters = -GE:not=8 -NW:check_nfs=no -OUTPUT:rtd=yes\n')
		ff.write('job = genome,denovo,accurate\n')
		if pair: 
			ff.write('readgroup = DataIlluminaPairedLib\n')
			ff.write('autopairing\n')
		else:
			ff.write('readgroup = DataIlluminaSingleLib\n')
		if pair: ff.write('data = '+wd+'/fastq/'+key+'_1mira.fastq '+ wd+'/fastq/'+key+'_2mira.fastq\n')
		else: ff.write('data = '+wd+'/fastq/'+key+'_1mira.fastq \n')
		ff.write('technology = solexa\n')
		ff.close()
		miraOut= wd+'/'+key+'_assembly/'+key+'_d_results/'+key+'_out.unpadded.fasta'
		rval.append(miraOut)
		job+=1
		if job%nservers==0: f.write('wait\n')
	f.close()
	return rval

def miraOLC(mira_local):
	f = open('miraOLC.sh', 'w')
	job=0
	rval=[]
	for (key, vdict) in seeds.items():
		key2=key+'OLC'
		serverTag = '{ time ssh '+servers[job%nservers]
		ff = open(wd+'/'+key2+'.conf', 'w')
		f.write(serverTag+' "cd '+wd+' && mira -t 8 '+wd+'/'+key2+'.conf" ; } 2> '+wd+'/'+key2+'_miraOLC.time &\n')
		ff.write('project = '+key2+'\n')
		if mira_local: ff.write('parameters = -GE:not=8 -DI:trt=/home/BSI/306307/ -OUTPUT:rtd=yes\n')
		else: ff.write('parameters = -GE:not=8 -NW:check_nfs=no -OUTPUT:rtd=yes \n')
		#ff.write('parameters = --noqualities=454\n')
		ff.write('job = genome,denovo,accurate\n')
		ff.write('readgroup\n')
		ff.write('technology = 454\n')
		ff.write('data = fa::'+wd+'/fastq/'+key+'_contig2 \n')
		#ff.write('data = '+wd+'/test.fasta \n')
		ff.write('default_qual=30\n')
		ff.close()
		miraOut= wd+'/'+key2+'_assembly/'+key2+'_d_results/'+key2+'_out.unpadded.fasta'
		rval.append(miraOut)
		job+=1
		if job%nservers==0: f.write('wait\n')
	f.close()
	return rval

def combineContig(contigLength1,r1,r2,r3,r4,r5,r6,r7):
#abyssOut= wd+'/'+key+'-unitigs.fa'
#velvetOut=wd+'/velvet_'+key+'/contigs.fa'
# /mnt/cluster/xdeng/script/statContigs.py /mnt/san2/xdeng/Vhunt/131210_TN_CDC5/soap_out/Phan491DNA_soap.contig soap Phan491DNA
# /mnt/cluster/xdeng/script/statContigs.py /mnt/san2/xdeng/Vhunt/131210_TN_CDC5/velvet_Phan491DNA/contigs.fa meta_velvet Phan491DNA
# /mnt/cluster/xdeng/script/statContigs.py /mnt/san2/xdeng/Vhunt/131210_TN_CDC5/abyss_Phan491DNA/abyss-unitigs.fa abyss Phan491DNA
# /mnt/cluster/xdeng/script/statContigs.py Phan491DNA_assembly/Phan491DNA_d_results/Phan491DNA_out.unpadded.fasta mira Phan491DNA
# /mnt/cluster/xdeng/script/faLenFilter.py /mnt/san2/xdeng/Vhunt/131210_TN_CDC5/fastq/Phan491DNA_contig /mnt/san2/xdeng/Vhunt/131210_TN_CDC5/fastq/Phan491DNA_contig2 150
	f = open(wd+'/combineContig.sh', 'w')
	i=0
	for (key, vdict) in seeds.items():
		combfile=r1[i]+' '+r2[i]+' '+r3[i]+' '+r4[i]+' '+r5[i]+' '+r6[i]+' '+r7[i]
		f.write('cat '+combfile+ ' > '+wd+'/fastq/'+key+'_contig\n')
		f.write(dirscr+'faLenFilter.py '+wd+'/fastq/'+key+'_contig '+wd+'/fastq/'+key+'_contig2 '+str(contigLength1) +'\n')
		i+=1
	f.close()

def combineContig2(contigLength1,assembly_para):
	f = open(wd+'/combineContig2.sh', 'w')
	combfile=[]
	contig_out = wd+'/contig_individual/'
	for (key, vdict) in seeds.items():
		if 'S' in assembly_para: combfile.append(contig_out+key+'_S.contig') #soap
		if 'V' in assembly_para: combfile.append(contig_out+key+'_V.contig') #velvet
		if 'A' in assembly_para: combfile.append(contig_out+key+'_A.contig') 
		if 'M' in assembly_para: combfile.append(contig_out+key+'_M.contig')
		if 'T' in assembly_para: combfile.append(contig_out+key+'_T.contig')
		if 's' in assembly_para: combfile.append(contig_out+key+'_s.contig') 
		if 'v' in assembly_para: combfile.append(contig_out+key+'_v.contig') 
		if 'a' in assembly_para: combfile.append(contig_out+key+'_a.contig') 
		f.write('cat '+' '.join(combfile)+ ' > '+wd+'/fastq/'+key+'_contig\n')
		combfile=[]
		f.write(dirscr+'faLenFilter.py '+wd+'/fastq/'+key+'_contig '+wd+'/fastq/'+key+'_contig2 '+str(contigLength1) +'\n')
	f.close()


def moveContig_Kmer(assembly_para, runid, r1,r2,r3,r4, r5,r6,r7):
	contig_out = wd+'/contig_'+runid+'/'
	try: os.mkdir(contig_out)
	except: pass
	f = open(wd+'/moveContig.sh', 'w')
	f1 = open(wd+'/statContig.sh', 'w')
	f2 = open(wd+'/blastContig.sh', 'w')
	f1.write('echo program label key top1 top2 top3 N50 ncontigs n500 n1000 n2000 n5000\n')
	i = 0
	for (key, vdict) in seeds.items():
		if 'S' in assembly_para: f.write('mv '+r1[i]+' '+contig_out+key+'_S.contig\n') #soap
		if 'V' in assembly_para:  f.write('mv '+r2[i]+' '+contig_out+key+'_V.contig\n') #velvet
		if 'A' in assembly_para: 
			f.write('cp --dereference '+r3[i]+' '+contig_out+key+'_A.contig\n')
			#f.write('rm '+r3[i]+'\n')#Abyss
		if 'M' in assembly_para: f.write('mv '+r4[i]+' '+contig_out+key+'_M.contig\n') #Mira
		if 's' in assembly_para: f.write('mv '+r5[i]+' '+contig_out+key+'_s.contig\n')
		if 'v' in assembly_para: f.write('mv '+r6[i]+' '+contig_out+key+'_v.contig\n')
		if 'a' in assembly_para: f.write('mv '+r7[i]+' '+contig_out+key+'_a.contig\n')
		if 'C' in assembly_para: f.write('mv '+wd+'/fastq/'+key+'_contig3 '+contig_out+key+'_C.contig\n')
		if 'O' in assembly_para: f.write('mv '+wd+'/fastq/'+key+'_contig2-contigs.fa '+contig_out+key+'_O.contig\n')
		if 'S' in assembly_para: f1.write(dirscr+'statContigs.py '+contig_out+key+'_S.contig  soap '+key+' '+runid+'\n')
		if 'V' in assembly_para: f1.write(dirscr+'statContigs.py '+contig_out+key+'_V.contig  meta '+key+' '+runid+'\n')
		if 'A' in assembly_para: f1.write(dirscr+'statContigs.py '+contig_out+key+'_A.contig  abyss '+key+' '+runid+'\n')
		if 'M' in assembly_para: f1.write(dirscr+'statContigs.py '+contig_out+key+'_M.contig  mira '+key+' '+runid+'\n')
		if 'C' in assembly_para: f1.write(dirscr+'statContigs.py '+contig_out+key+'_C.contig cap3 '+key+' '+runid+'\n')
		if 'O' in assembly_para: f1.write(dirscr+'statContigs.py '+contig_out+key+'_O.contig minimo '+key+' '+runid+'\n')
		if 'S' in assembly_para: f2.write(dirscr+'blastContig.py '+contig_out+key+'_S.contig  soap_'+soapkmer+'_'+key+' '+blastRef+'\n') #use BASV for ref 5 referneces
		if 'V' in assembly_para: f2.write(dirscr+'blastContig.py '+contig_out+key+'_V.contig  meta_'+metakmer+'_'+key+' '+blastRef+'\n') # or BASV
		if 'A' in assembly_para: f2.write(dirscr+'blastContig.py '+contig_out+key+'_A.contig  abyss_'+abysskmer+'_'+key+' '+blastRef+'\n')
		if 'M' in assembly_para: f2.write(dirscr+'blastContig.py '+contig_out+key+'_M.contig  mira_'+key+' '+blastRef+'\n')
		if 'C' in assembly_para: f2.write(dirscr+'blastContig.py '+contig_out+key+'_C.contig cap3_'+key+' '+blastRef+'\n')
		if 'O' in assembly_para: f2.write(dirscr+'blastContig.py '+contig_out+key+'_O.contig minimo_'+key+' '+blastRef+'\n')
		i+=1
	f1.close()
	f2.close()
	# f.write('mv '+wd+'/assembly.sh '+contig_out+'\n')
	f.write('mv '+wd+'/*.time '+contig_out+'\n') 
	f.close()

#move contig for individual
def moveContig1(assembly_para, r1,r2,r3,r4, r5,r6,r7, r8, r9, r10, r11, r12):
	contig_out = wd+'/contig_individual/'
	try: os.mkdir(contig_out)
	except: pass
	f = open(wd+'/moveContig.sh', 'w')
	f1 = open(wd+'/statContig.sh', 'w')
	f2 = open(wd+'/blastContig.sh', 'w')
	f1.write('echo program label key top1 top2 top3 N50 ncontigs n500 n1000 n2000 n5000\n')
	i = 0
	for (key, vdict) in seeds.items():
		if 'S' in assembly_para: f.write('mv '+r1[i]+' '+contig_out+key+'_S.contig\n') #soap
		if 'V' in assembly_para:  f.write('mv '+r2[i]+' '+contig_out+key+'_V.contig\n') #velvet
		if 'A' in assembly_para: 
			f.write('cp --dereference '+r3[i]+' '+contig_out+key+'_A.contig\n')
			#f.write('rm '+r3[i]+'\n')#Abyss
		if 'M' in assembly_para: f.write('mv '+r4[i]+' '+contig_out+key+'_M.contig\n') #Mira
		if 's' in assembly_para: f.write('mv '+r5[i]+' '+contig_out+key+'_s.contig\n')
		if 'v' in assembly_para: f.write('mv '+r6[i]+' '+contig_out+key+'_v.contig\n')
		if 'a' in assembly_para: f.write('mv '+r7[i]+' '+contig_out+key+'_a.contig\n')
		if 'G' in assembly_para: f.write('mv '+r8[i]+' '+contig_out+key+'_G.contig\n') # omega
		if 'W' in assembly_para: f.write('mv '+r9[i]+' '+contig_out+key+'_W.contig\n')  #celera
		if 'T' in assembly_para: f.write('mv '+r10[i]+' '+contig_out+key+'_T.contig\n')  #celera
		if 'X' in assembly_para: f.write('mv '+r11[i]+' '+contig_out+key+'_X.contig\n')  #celera
		if 'I' in assembly_para: f.write('mv '+r12[i]+' '+contig_out+key+'_I.contig\n')  #idba
		
		if 'S' in assembly_para: f1.write(dirscr+'statContigs.py '+contig_out+key+'_S.contig  soap '+key+' '+assembly_para+'\n')
		if 'V' in assembly_para: f1.write(dirscr+'statContigs.py '+contig_out+key+'_V.contig  meta '+key+' '+assembly_para+'\n')
		if 'A' in assembly_para: f1.write(dirscr+'statContigs.py '+contig_out+key+'_A.contig  abyss '+key+' '+assembly_para+'\n')
		if 'M' in assembly_para: f1.write(dirscr+'statContigs.py '+contig_out+key+'_M.contig  mira '+key+' '+assembly_para+'\n')
		if 'G' in assembly_para: f1.write(dirscr+'statContigs.py '+contig_out+key+'_G.contig  omega '+key+' '+assembly_para+'\n')
		if 'W' in assembly_para: f1.write(dirscr+'statContigs.py '+contig_out+key+'_W.contig  celera '+key+' '+assembly_para+'\n')
		if 'T' in assembly_para: f1.write(dirscr+'statContigs.py '+contig_out+key+'_T.contig  trinity '+key+' '+assembly_para+'\n')
		if 'X' in assembly_para: f1.write(dirscr+'statContigs.py '+contig_out+key+'_X.contig  masurca '+key+' '+assembly_para+'\n')
		if 'I' in assembly_para: f1.write(dirscr+'statContigs.py '+contig_out+key+'_I.contig  idba '+key+' '+assembly_para+'\n')
		
		if 'S' in assembly_para: f2.write(dirscr+'blastContig.py '+contig_out+key+'_S.contig  soap_'+soapkmer+'_'+key+' '+blastRef+'\n')
		if 'V' in assembly_para: f2.write(dirscr+'blastContig.py '+contig_out+key+'_V.contig  meta_'+metakmer+'_'+key+' '+blastRef+'\n')
		if 'A' in assembly_para: f2.write(dirscr+'blastContig.py '+contig_out+key+'_A.contig  abyss_'+abysskmer+'_'+key+' '+blastRef+'\n')
		if 'M' in assembly_para: f2.write(dirscr+'blastContig.py '+contig_out+key+'_M.contig  mira_'+key+' '+blastRef+'\n')
		if 'G' in assembly_para: f2.write(dirscr+'blastContig.py '+contig_out+key+'_G.contig  omega_'+key+' '+blastRef+'\n')
		if 'W' in assembly_para: f2.write(dirscr+'blastContig.py '+contig_out+key+'_W.contig  celera_'+key+' '+blastRef+'\n')
		if 'T' in assembly_para: f2.write(dirscr+'blastContig.py '+contig_out+key+'_T.contig  trinity_'+key+' '+blastRef+'\n')
		if 'X' in assembly_para: f2.write(dirscr+'blastContig.py '+contig_out+key+'_X.contig  masurca_'+key+' '+blastRef+'\n')
		if 'I' in assembly_para: f2.write(dirscr+'blastContig.py '+contig_out+key+'_I.contig  idba_'+key+' '+blastRef+'\n')
		
		i+=1
	f1.close()
	f2.close()
	# f.write('mv '+wd+'/assembly.sh '+contig_out+'\n')
	f.write('mv '+wd+'/*.time '+contig_out+'\n') 
	f.close()

#move contig for ensembled assemblers
def moveContig2(assembly_para):
	contig_out = wd+'/contig_'+assembly_para+'/'
	try: os.mkdir(contig_out)
	except: pass
	f = open(contig_out+'moveContig.sh', 'w')
	f1 = open(contig_out+'statContig.sh', 'w')
	f2 = open(contig_out+'blastContig.sh', 'w')
	f1.write('echo program label key top1 top2 top3 N50 ncontigs n500 n1000 n2000 n5000\n')
	i = 0
	for (key, vdict) in seeds.items():
		key2=key+'OLC'
		trinityOut = wd+'/trinity_'+key+'trinityOLC'+'/'+'Trinity.fasta'
		if 'C' in assembly_para: f.write('mv '+wd+'/fastq/'+key+'_contig3 '+contig_out+key+'_C.contig\n')
		if 'O' in assembly_para: f.write('mv '+wd+'/fastq/'+key+'_contig2-contigs.fa '+contig_out+key+'_O.contig\n')
		if 'R' in assembly_para: f.write('mv '+wd+'/'+key2+'_assembly/'+key2+'_d_results/'+key2+'_out.unpadded.fasta '+contig_out+key+'_R.contig\n')
		if 'Y' in assembly_para: f.write('mv '+trinityOut+' '+contig_out+key+'_Y.contig\n')
		if 'C' in assembly_para: f1.write(dirscr+'statContigs.py '+contig_out+key+'_C.contig cap3 '+key+' '+assembly_para+'\n')
		if 'O' in assembly_para: f1.write(dirscr+'statContigs.py '+contig_out+key+'_O.contig minimo '+key+' '+assembly_para+'\n')
		if 'R' in assembly_para: f1.write(dirscr+'statContigs.py '+contig_out+key+'_R.contig miraOLC '+key+' '+assembly_para+'\n')
		if 'Y' in assembly_para: f1.write(dirscr+'statContigs.py '+contig_out+key+'_Y.contig trinityOLC '+key+' '+assembly_para+'\n')
		if 'C' in assembly_para: f2.write(dirscr+'blastContig.py '+contig_out+key+'_C.contig cap3_'+key+' '+blastRef+'\n')
		if 'O' in assembly_para: f2.write(dirscr+'blastContig.py '+contig_out+key+'_O.contig minimo_'+key+' '+blastRef+'\n')
		if 'R' in assembly_para: f2.write(dirscr+'blastContig.py '+contig_out+key+'_R.contig miraOLC_'+key+' '+blastRef+'\n')
		if 'Y' in assembly_para: f2.write(dirscr+'blastContig.py '+contig_out+key+'_Y.contig trinityOLC_'+key+' '+blastRef+'\n')
		i+=1
	f1.close()
	f2.close()
	# f.write('mv '+wd+'/assembly.sh '+contig_out+'\n')
	f.write('mv '+wd+'/*.time '+contig_out+'\n') 
	f.write('rm -rf '+wd+'/fastq/*contig* \nwait\n')
	f.close()

# def cap3():
# #ssh bsidna3 cap3 /mnt/san2/xdeng/Vhunt/131210_TN_CDC5/fastq/Phan491DNA_contig2 &
	# f = open('cap3.sh', 'w')
	# job=0
	# for (key, vdict) in seeds.items():
		# serverTag = '{ time ssh '+servers[job%nservers]
		# f.write(serverTag+' "cap3 '+wd+'/fastq/'+key+'_contig2" ; } 2> '+wd+'/'+key+'_cap3.time &\n')
		# job+=1
		# if job%nservers==0: f.write('wait\n')
	# f.close()

def cap3():
#ssh bsidna3 cap3 /mnt/san2/xdeng/Vhunt/131210_TN_CDC5/fastq/Phan491DNA_contig2 &
	f = open('cap3.sh', 'w')
	job=0
	for (key, vdict) in seeds.items():
		serverTag = '{ time ssh '+servers[job%nservers]
		f.write(serverTag+' "cap3 '+wd+'/fastq/'+key+'_contig2" ; } 2> '+wd+'/'+key+'_cap3.time &\n')
		job+=1
		if job%nservers==0: f.write('wait\n')
	f.close()

def minimo():
#ssh bsidna3 cap3 /mnt/san2/xdeng/Vhunt/131210_TN_CDC5/fastq/Phan491DNA_contig2 &
	f = open('minimo.sh', 'w')
	job=0
	for (key, vdict) in seeds.items():
		serverTag = '{ time ssh '+servers[job%nservers]
		f.write(serverTag+' "'+Minimo+' '+wd+'/fastq/'+key+'_contig2 -D FASTA_EXP=1" ; } 2> '+wd+'/'+key+'_minimo.time &\n')
		job+=1
		if job%nservers==0: f.write('wait\n')
	f.close()

# def minimoFilter(contigLength2):
	# f = open('minimoFilter.sh', 'w')
	# job=0
	# for (key, vdict) in seeds.items():
		# #serverTag = 'ssh '+servers[job%nservers]
		# serverTag=''
		# f.write(serverTag+' '+dirscr+'faLenFilter.py '+wd+'/fastq/'+key+'_contig2-contigs.fa '+wd+'/fastq/'+key+'_contig4 '+str(contigLength2)+' '+base+'_'+key+'\n')
		# # job+=1
		# # if job%nservers==0: f.write('wait\n')
	# f.close()

def trinityOLC(RAM='25G'):
	f = open('trinityOLC.sh', 'w')
	job=0
	for (key, vdict) in seeds.items():
		key1 = key+'trinityOLC'
		outdir = wd+'/trinity_'+key1+'/'
		trinityOut=wd+'/trinity_'+key+'trinityOLC'+'/'+'Trinity.fasta'
		serverTag = '{ time ssh '+servers[job%nservers]
		file1 = wd+'/fastq/'+key+'_contig2'
		f.write(serverTag+' "'+trinitypath+'/Trinity --seqType fa --JM '+RAM+' --output '+outdir+' --single '+file1 + \
			' --CPU 8 "  ; } 2> '+wd+'/'+key1+'_trinityOLC.time &\n')
		job+=1
		if job%nservers==0: f.write('wait\n')
	f.close()
	
# def TrinityOLCFilter(contigLength2):
	# f = open('minimoFilter.sh', 'w')
	# job=0
	# for (key, vdict) in seeds.items():
		# #serverTag = 'ssh '+servers[job%nservers]
		# serverTag=''
		# key1 = key+'OLC'
		# f.write(serverTag+' '+dirscr+'faLenFilter.py '+wd+'/trinity_'+key1+'/Trinity.fasta '+wd+'/fastq/'+key+'_contig4 '+str(contigLength2)+' '+base+'_'+key+'\n')
		# # job+=1
		# # if job%nservers==0: f.write('wait\n')
	# f.close()

# def miraOLCFilter(contigLength2):
	# f = open('miraOLCFilter.sh', 'w')
	# job=0
	# for (key, vdict) in seeds.items():
		# #serverTag = 'ssh '+servers[job%nservers]
		# serverTag=''
		# key2=key+'OLC'
		# f.write(serverTag+' '+dirscr+'faLenFilter.py '+wd+'/'+key2+'_assembly/'+key2+'_d_results/'+key2+'_out.unpadded.fasta '+wd+'/fastq/'+key+'_contig4 '+str(contigLength2)+' '+base+'_'+key+'\n')
		# # job+=1
		# # if job%nservers==0: f.write('wait\n')
	# f.close()

def cap3MergeSinglet(contigLength2):
	f = open('cap3MergeSinglet.sh', 'w')
	job=0
	for (key, vdict) in seeds.items():
		#serverTag = 'ssh '+servers[job%nservers]
		serverTag=''
		f.write(serverTag+' cat '+wd+'/fastq/'+key+'_contig2.cap.singlets '+wd+'/fastq/'+key+'_contig2.cap.contigs > '+ wd+'/fastq/'+key+'_contig3\n')
		f.write(serverTag+' '+dirscr+'faLenFilter.py '+wd+'/fastq/'+key+'_contig3 '+ wd+'/fastq/'+key+'_contig4 '+str(contigLength2)+' '+base+'_'+key+'\n')
		# job+=1
		# if job%nservers==0: f.write('wait\n')
	f.close()

def Test_Kmer(assembly_para, abysskmer, soapkmer, metakmer):
	runid = assembly_para+'_'+soapkmer+'_'+abysskmer+'_'+metakmer
	contig_out = wd+'/contig_'+runid+'/'
	try: os.mkdir(contig_out)
	except: pass
	sf = open('assembly'+soapkmer+'.sh', 'w')
	sf.write('echo Welcome13 |sudo -S chmod 777 * -R\n')
	sf.write('rm -rf soap_out/* \nwait\n')
	sf.write('rm -rf *.time \nwait\n')
	sf.write('rm -rf velvet_* \nwait\n')
	sf.write('rm -rf abyss_*/* \nwait\n')
	if pair: r1=soap_pair()
	else: r1=soap_single()
	sf.write('source '+contig_out+'soap.sh\nwait\n')
	r2=meta_velvet()
	sf.write('source '+contig_out+'meta_velvet.sh\nwait\n')
	r3=abyss()
	sf.write('source '+contig_out+'abyss.sh\nwait\n')
	os.system('mv soap.sh '+contig_out)
	os.system('mv abyss.sh '+contig_out)
	os.system('mv meta_velvet.sh '+contig_out)

	r4=mira4(mira_local)
	r5=soap_partition()
	r6=velvet_partition()
	r7=abyss_partition()

	moveContig_Kmer(assembly_para, runid, r1,r2,r3,r4,r5,r6,r7)
	os.system('mv moveContig.sh '+contig_out)
	os.system('mv statContig.sh '+contig_out)
	os.system('mv blastContig.sh '+contig_out)
	sf.write('echo Welcome13 |sudo -S chmod 777 * -R\n')
	sf.write('source '+contig_out+'moveContig.sh\n')
	sf.write('source '+contig_out+'statContig.sh > '+contig_out+'contig.log\n')
	sf.write('source blastContig.sh > '+contig_out+'contig.blast\n')
	sf.close()

def Assembly_individual(assembly_para):
	sf = open('assembly_individual.sh', 'w')
	partition()
	if ('s' in assembly_para) or ('v' in assembly_para) or ('a' in assembly_para):
		sf.write('source partition.sh\nwait\n')
	if 'S' in assembly_para: 
		if pair: r1=soap_pair()
		else: r1=soap_single()
		sf.write('source soap.sh\nwait\n')
	if 'V' in assembly_para: 
		r2=meta_velvet()
		sf.write('source meta_velvet.sh\nwait\n')
	if 'A' in assembly_para: 
		r3=abyss()
		sf.write('source abyss.sh\nwait\n')
	if 'M' in assembly_para: 
		r4=mira4(mira_local)
		sf.write('source mira.sh\nwait\n')
	if 'G' in assembly_para: 
		r8=omega()
		sf.write('source omega.sh\nwait\n')
	if 'T' in assembly_para: 
		r10=trinity(RAM)
		sf.write('source trinity.sh\nwait\n')
	if 'X' in assembly_para: 
		r11=masurca()
		sf.write('source masurca.sh\nwait\n')
	if 'I' in assembly_para: 
		r12=idba()
		sf.write('source idba.sh\nwait\n')
	if 'W' in assembly_para: 
		r9=celera()
		sf.write('source celera_prepReads.sh\nwait\n')
		sf.write('source celeraPrep.sh\nwait\n')
		sf.write('source celera.sh\nwait\n')
		sf.write('source celeraCombine.sh\nwait\n')
	if 's' in assembly_para:
		r5=soap_partition()
		sf.write('source soap_partition.sh\nwait\n')
		sf.write('source soap_combine.sh\nwait\n')
	if 'v' in assembly_para:
		r6=velvet_partition()
		sf.write('source velvet_partition.sh\nwait\n')
		sf.write('source velvet_combine.sh\nwait\n')
	if 'a' in assembly_para:
		r7=abyss_partition()
		sf.write('source abyss_partition.sh\nwait\n')
		sf.write('source abyss_combine.sh\nwait\n')
	moveContig1(assembly_para, r1,r2,r3,r4, r5,r6,r7,r8,r9,r10,r11, r12)
	sf.write('echo Welcome13 |sudo -S chmod 777 * -R\n')
	sf.write('source moveContig.sh\n')
	contig_out = wd+'/contig_individual/'
	sf.write('source statContig.sh > '+contig_out+'contig.log\n')
	sf.write('source blastContig.sh > '+contig_out+'contig.blast\n')
	sf.write('source abyss_partition_blast.sh > abyss_chunk_contigblast.chunk\n')
	sf.close()

def Assembly_combine(assembly_para):
	contig_out = wd+'/contig_'+assembly_para+'/'
	try: os.mkdir(contig_out)
	except: pass

	sf = open(contig_out+'assembly_combine.sh', 'w')
	combineContig2(contigLength1,assembly_para)
	os.system('mv combineContig2.sh '+contig_out+'\nwait\n')
	sf.write('source '+contig_out+'combineContig2.sh \nwait\n')
	if 'C' in assembly_para: 
		cap3()
		cap3MergeSinglet(contigLength2)
		os.system('mv cap3.sh '+contig_out+'\nwait\n')
		os.system('mv cap3MergeSinglet.sh '+contig_out+'\nwait\n')
		sf.write('source '+contig_out+'cap3.sh\nwait\n')
		sf.write('source '+contig_out+'cap3MergeSinglet.sh \nwait\n')
	elif 'O' in assembly_para:
		minimo()
		#minimoFilter(contigLength2)
		os.system('mv minimo.sh '+contig_out+'\nwait\n')
		#os.system('mv minimoFilter.sh '+contig_out+'\nwait\n')
		sf.write('source '+contig_out+'minimo.sh\nwait\n')
		#sf.write('source '+contig_out+'minimoFilter.sh\nwait\n')
	elif 'R' in assembly_para:
		miraOLC(mira_local)
		# miraOLCFilter(contigLength2)
		os.system('mv miraOLC.sh '+contig_out+'\nwait\n')
		#os.system('mv miraOLCFilter.sh '+contig_out+'\nwait\n')
		sf.write('source '+contig_out+'miraOLC.sh\nwait\n')
		#sf.write('source '+contig_out+'miraOLCFilter.sh\nwait\n')
	elif 'Y' in assembly_para:
		trinityOLC(RAM)
		os.system('mv trinityOLC.sh '+contig_out+'\nwait\n')
		sf.write('source '+contig_out+'trinityOLC.sh\nwait\n')

	moveContig2(assembly_para)
	sf.write('echo Welcome13 |sudo -S chmod 777 * -R\n')
	sf.write('source '+contig_out+'moveContig.sh\n')
	sf.write('source '+contig_out+'statContig.sh > '+contig_out+'contig.log\n')
	sf.write('source '+contig_out+'blastContig.sh > '+contig_out+'contig.blast\n')
	sf.close()

def AddPipe(para, sf):
	assembly_para=para
	Assembly_combine(assembly_para)
	sf.write('source '+wd+'/contig_'+assembly_para+'/'+'assembly_combine.sh\n')

if __name__ == "__main__":
	blastRef = sys.argv[1] #'NIBSC' #'ref5Virus' #'BASV'# 'SURPI' #'bacteria' #'BASV' #'ref5' # BASV, NIBSC, bacteria
	try: RAM = sys.argv[2]
	except: RAM='25G'
	sra=False
	skipall=True #don't do any cleaning
	keep_human=True
	mira_local=True
	keep_bac=True
	pair=True
	skipadaptor=True
	skipread=False
	thread = '8'
	abysskmer='31'
	soapkmer='31'
	metakmer='31'
	hsp='NO' # show whole reads, if 'YES' show only blast hsp
	wd = os.path.abspath(os.path.dirname('.'))
	#genseedfile()
	if pair:
		seeds=readSeeds2()
	else:
		seeds=readSeeds1()
	if not wd.startswith('/mnt/'): wd='/mnt'+wd
	base = os.path.basename(wd)
	print 'blastRef', blastRef
	print 'path',wd
	print 'thread', thread
	print 'base', base
	for key in seeds.keys():
		print key, seeds[key]
	print 'keep_human', keep_human
	print 'keep_bac', keep_bac
	print 'pairend', pair
	print 'skipall', skipall
	print 'skipadaptor', skipadaptor
	print 'skipread', skipread
	print 'hsp', hsp
	print 'metakmer', metakmer
	print 'soapkmer', soapkmer
	print 'abysskmer', abysskmer
	print 'servers', servers
	oof=open(wd+'/run.log', 'w')
	for key in seeds.keys():
		print >> oof, key, seeds[key]
	print >>oof, 'keep_human', keep_human
	print >>oof, 'keep_bac', keep_bac
	print >>oof,'pairend', pair
	print >>oof,'skipadaptor', skipadaptor
	print >>oof,'skipread', skipread
	print >>oof,'hsp', hsp
	print >>oof, 'metakmer', metakmer
	print >>oof, 'soapkmer', soapkmer
	print >>oof, 'abysskmer', abysskmer
	oof.close()
	# sampleFastq()
	# sys.exit()
	n=140 #number fo splits
	length=49 #lengh reads for blast
	contigLength1 = 300 #length before cap3
	contigLength2 = 300 #length after cap3 before blast
	myslen = 1200
	print 'read threshold', length
	print 'contig before cap3', contigLength1
	print 'contig after cap3', contigLength2
	print 'mys threshold', myslen
	prep_sra()
	polyA()
	# bowtieHuman(pair)
	# #bowtieNT()
	skipbowtie()
	bowtieBac(pair, keep_human)
	sam2fq(pair, keep_human)

	check_pair_overlap()
	prepBlastFile_adaptor()
	trim()
	skip_adaptor()
	fq_check()
	prep_reads(n, length, pair)
	clean_dir()
	assembly_para='SAV' #SAOP, Abyss, Velvet, Mira, Cap3, no Minimo
	#assembly_para='SAVMsavmO' #SAOP, Abyss, Velvet, Mira, Cap3, no Minimo
	bowtieRef(blastRef)
	# abysskmer='31'; soapkmer='31'; metakmer='31'
	# Test_Kmer(assembly_para, abysskmer, soapkmer, metakmer)
	# abysskmer='41'; soapkmer='41'; metakmer='41'
	# Test_Kmer(assembly_para, abysskmer, soapkmer, metakmer)
	# abysskmer='51'; soapkmer='51'; metakmer='51'
	# Test_Kmer(assembly_para, abysskmer, soapkmer, metakmer)
	# abysskmer='61'; soapkmer='61'; metakmer='61'
	# Test_Kmer(assembly_para, abysskmer, soapkmer, metakmer)
	
	#if pair: prepPriceFile() #not doing price
	sf=open('pipeline_run.sh', 'w')
	# os.system('PriceTI  -fp s1.fq s2.fq 300 -icf seeds.fa 1 1 5  -nc 30 -dbmax 72 -mol 35 -tol 20 -mpi 80 -target 90 2 2 2 -a 7 -o price.fa')

	if sra: 
		sf.write('source sra.sh\n')
	if not skipall:
		if keep_human and keep_bac:
			sf.write('source skipbowtie.sh\n')
		else:
			sf.write('source bowtieBac.sh\nwait\n')
			sf.write('source bowtiesam2fq.sh  >>stats.log\nwait\n')

		sf.write('source clonetrim.sh >>stats.log \nwait\n')
		if skipadaptor:
			sf.write('source skipadaptor.sh\n')
		else:
			sf.write('source clonefq2fa.sh \nwait\n')
			sf.write('source prepBlastFile_adaptor.sh \n')
			sf.write('source blast_adaptor.sh \nwait\n')
			sf.write('source blasttrim.sh >>stats.log \nwait\n')
			sf.write('source qualitytrim.sh >>stats.log \nwait\n')
		sf.write('source fq_clean.sh \nwait\n')
	else:
		pass
	#sf.write('source check_pair.sh >>stats.log \nwait\n')
	# sf.write('source polyA_raw.sh >>stats.log \n')
	# sf.write('source polyA_clean.sh >>stats.log \n')
	sf.write('source prep_reads.sh\nwait\n')
	
	# #kmer test
	# sf.write('source assembly31.sh\nwait\n')
	# sf.write('source assembly41.sh\nwait\n')
	# sf.write('source assembly51.sh\nwait\n')
	# sf.write('source assembly61.sh\nwait\n')
	# sf.close()
	# sys.exit()
	# sf.write('source clean.sh\nwait\n')
	abysskmer='31'
	soapkmer='31'
	metakmer='31'
	assembly_para='SAVMGWTXIsav' #SAOP, Abyss, Velvet, Mira, Cap3, no Minimo
	Assembly_individual(assembly_para)
	sf.write('source assembly_individual.sh\n')
	
	#assembly emsemble
	AddPipe('SC', sf)
	AddPipe('AC', sf)
	AddPipe('VC', sf)
	AddPipe('SO', sf)
	AddPipe('AO', sf)
	AddPipe('VO', sf)
	AddPipe('SsC', sf)
	AddPipe('AaC', sf)
	AddPipe('VvC', sf)
	AddPipe('SsO', sf)
	AddPipe('AaO', sf)
	AddPipe('VvO', sf)
	
	AddPipe('SAVC', sf)
	AddPipe('SAVMC', sf)
	AddPipe('SAVO', sf)
	AddPipe('SAVMO', sf)
	
	AddPipe('SAVaC', sf)
	AddPipe('SAVaO', sf)
	AddPipe('SAVMaO', sf)
	AddPipe('SAVMaC', sf)
	AddPipe('SAVaR', sf)
	AddPipe('SAVTaC', sf)
	AddPipe('SAVaY', sf) #trinity
	#bowtieBASV()
	sf.write('grep ""  contig_*/contig.blast > contig.blast.all\n')
	# sf.write('tail -n +1  contig_*/contig.blast > contig.blast.all.html\n')
	sf.write('contightml.py contig.blast.all contig.blast.all.html\n')
	sf.write('tail -n 3 contig_*/*.time > contig.timing\n')
	sf.write('timing.py contig.timing > contig.timing2\n')
	sf.write('cat contig*/contig.log > contig.log\n')
	sf.write('svg_blast.py contig.blast.all > contig.blast.table\n')
	# sf.write('source bowtieRef.sh\n')
	sf.close()
