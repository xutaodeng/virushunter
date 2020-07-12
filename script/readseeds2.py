#!/usr/bin/env python
from collections import defaultdict
import os.path
import os
import string
import sys

servers=['bsidna3','bsidna4','bsidna5','bsidna6','bsidna7','bsidna8','bsidna9',\
'bsidna11','bsidna12','bsidna13','bsidna14','bsidna15','bsidna16','bsidna17','bsidna18','bsidna19','bsidna20',\
'bsidna21','bsidna22','bsidna23','bsidna24','bsidna25','bsidna27','bsidna29', 'bsidna30']
soappath='/mnt/cluster/tools/SOAPdenovo-63mer'
velvetg='/mnt/cluster/tools/velvet_1.2.10/velvetg'
velveth='/mnt/cluster/tools/velvet_1.2.10/velveth'
meta_velvetg="/mnt/cluster/tools/MetaVelvet-1.2.02/meta-velvetg"
bowtiepath='/mnt/san/cluster/tools/bowtie2-2.1.0/bowtie2' 
blastxpath='/mnt/cluster/tools/ncbi-blast-2.2.27+/bin/blastx'
bowtieindexpath='/mnt/cluster/xdeng/hg19/mrnadna_bowtie'
bowtiebacs= ['/mnt/san/cluster/xdeng/nt/Bacteria1', \
			'/mnt/san/cluster/xdeng/nt/Bacteria2', \
			'/mnt/san/cluster/xdeng/nt/Bacteria3', \
			'/mnt/san/cluster/xdeng/nt/Bacteria4']
virusdbpath='/mnt/cluster/xdeng/blastdb/virus/virus_mask'
nvnrdbpath='/mnt/san/cluster/xdeng/blastdb/refseq/nvrefseq'
dirscr='/mnt/cluster/xdeng/script/'
Raytool='/mnt/cluster/tools/Ray2.3/Ray'
picard = '/mnt/cluster/tools/picard-tools-1.105/'
nservers=len(servers)

def genseedfile():
	os.system('cd fastq && ls -1 *fastq.gz > samples.txt && cd ..  ')

# def readSeeds1(seedfile='fastq/samples.txt'):
	# seeds=defaultdict(list)
	# f=open(seedfile, 'r')
	# for line in f:
		# key2=line.strip().rsplit('.', 1)[0]
		# key1=line.strip().split('_')[0]
		# if len(key1)<len(key2):
			# key=key1
		# else:
			# key=key2
		# seeds[key].append(line.strip())
	# f.close()
	# return seeds

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

def fa2fq():
	f = open('fa2fq.sh', 'w')
	for (key, fqfiles) in seeds.items():
		for fqfile in fqfiles:
			fafile = fqfile.replace('fastq', 'fasta')
			f.write('fa2fq2.py '+wd+'/fastq/'+fafile+' '+wd+'/fastq/'+fqfile+'\n') #quality trimming
	f.close()
	
def bam2fq():
	f = open('bam2fq.sh', 'w')
	for (key, fqfiles) in seeds.items():
		for fqfile in fqfiles:
			bamfile = fqfile.replace('fastq', 'bam')
			sam2fq = 'java  -Xms1024m -Xmx8g -jar  '+picard+'SamToFastq.jar INPUT='+wd+'/fastq/'+bamfile+' FASTQ='+wd+'/fastq/'+fqfile+' VALIDATION_STRINGENCY=SILENT'
			f.write(sam2fq+'\n')
	f.close()
	

def skip_adaptor():
	f = open('skipadaptor.sh', 'w')
	for (key, fqfiles) in seeds.items():
		i=1
		for fqfile in fqfiles:
		#ff = open(wd+'/soap_config/'+key+'_soap.config', 'w')
			index=key+'_'+str(i)
			fid = wd+'/fastq/'+index
			f.write('recodeID.py '+fid+'.dup '+' '+wd+'/fastq/' + key+'_'+str(i)+'_sequence.txt '+base+'_'+key+ ' '+str(i)+'\n') #quality trimming
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
			f.write(serverTag + ' '+dirscr+'clone_reads_rm_pair.py '+fqfil+' '+fid+'.dup &\n')
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
	f2.write('find ./fastq/ -type f ! -name "*.gz"  -a ! -name "*.fastq" ! -name "*.log" ! -name "*.php" ! -name "*.gif" -delete\nwait\n')
	f2.write('rm -rf blast_* \nwait\n')
	f2.write('rm -rf soap* \nwait\n')
	f2.write('rm -rf run \nwait\n')
	f2.write('rm -rf *.sh \nwait\n')
	f2.write('rm -rf *.log \nwait\n')
	f2.close()

def polyA():
	job=0
	f = open('polyA_raw.sh', 'w')
	f2 = open('polyA_clean.sh', 'w')
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

def prepBlastFile():
	try: os.mkdir(wd+'/'+base+'/blast/')
	except: pass
	f = open('prepBlastFile.sh', 'w')
	path='fastq/'
	f.write('cd fastq/ \n')
	bpath=wd+'/'+base+'/blast/'
	directory = os.path.basename(wd)
	alldb=[]
	for (key, fqfiles) in seeds.items():
		f.write('makeblastdb -dbtype nucl -parse_seqids -in '+key+'.fa -out '+key+'_blastdb\n')
		alldb.append(key+'_blastdb')
	f.write('blastdb_aliastool -dblist '+'\"'+' '.join(alldb)+'\" -dbtype nucl -out '+directory+'_blastdb -title \"'+directory+'_blastdb\"\n')
	f.write('mv *_blastdb* '+bpath+'\n')
	f.write('cd ..\n')
	f.close()

def prepPriceFile():
	try: os.mkdir(wd+'/'+base+'/price/')
	except: pass
	path=wd+'/'+base+'/price/'
	f = open('prepPriceFile.sh', 'w')
	directory = os.path.basename(wd)
	for (key, fqfiles) in seeds.items():
		i=1
		infq=[]
		outfq=[]
		for fqfile in fqfiles:
			index=key+'_'+str(i)
			outfq.append(path+index+'.fq')
			infq.append('fastq/'+index +'_sequence.txt ')
			i+=1
		f.write(dirscr+'fqprice.py '+infq[0]+' '+infq[1]+' '+outfq[0]+' '+outfq[1]+'\n')
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

# def bowtieHuman(pair):
	# f = open('bowtieHuman.sh', 'w')
	# f2 = open('bowtiesam2fq.sh', 'w')
	# job=0
	# job2=0
	# for (key, fqfiles) in seeds.items():
		# i=1
		# sams=[]
		# fqfils=[]
		# for fqfile in fqfiles:
			# fastqfile=wd+'/fastq/'+fqfile
			# index=key+'_'+str(i)
			# sam=wd+'/fastq/'+ index+'_0.sam'
			# sams.append(sam)
			# fqfil=wd+'/fastq/'+index+'.fil'
			# fqfils.append(fqfil)
			# serverTag = 'ssh '+servers[job%nservers]
			# f.write(serverTag+' '+bowtiepath+' --quiet --local --no-hd --reorder -p 7 -x '+bowtieindexpath+' -U '+fastqfile+' -S '+sam+' &\n')
			# job+=1
			# if job%nservers==0: f.write('wait\n')
			# # f.write('/mnt/cluster/tools/bwa-0.5.9/bwa aln -n 0.2 -o 3 -t 7 /mnt/cluster/xdeng/hg19/mrnadna.fa fastq/' +index +'_sequence.txt.tmp > fastq/'+ index+'.sai \n')
			# # f.write('/mnt/cluster/tools/bwa-0.5.9/bwa samse /mnt/cluster/xdeng/hg19/mrnadna.fa fastq/' + index+'.sai fastq/'+ index+'_sequence.txt.tmp > fastq/'+ index+'.sam \n')
			# i+=1
		# serverTag2 = 'ssh '+servers[job2%nservers]
		# if pair: f2.write(serverTag2 + ' '+dirscr+'sam2fq.py '+sams[0]+' '+sams[1]+' '+fqfils[0]+' '+fqfils[1]+ ' & \n')
		# else: f2.write(serverTag2 + ' '+dirscr+'sam2fq.py '+sams[0]+' '+fqfils[0]+ ' & \n')
		# job2+=1
		# if job2%nservers==0: f2.write('wait\n')
	# f.close()
	# f2.close()

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
			f.write(meta_velvetg+' '+outdir+' -ins_length 260"  ; } 2> '+key+'_meta.time &\n')
		else:
			f.write(serverTag+' "'+velveth+' '+outdir+' '+metakmer)
			f.write(' -short -fastq '+wd+'/fastq/'+key+'_1_sequence.txt && ')
			f.write(velvetg+' '+outdir+' -exp_cov auto && ')
			f.write(meta_velvetg+' '+outdir+'"  ; } 2> '+key+'_soap.time &\n')
		velvetOut=wd+'/velvet_'+key+'/contigs.fa'
		rval.append(velvetOut)
		job+=1
		if job%nservers==0: f.write('wait\n')
	f.close()
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
		f.write(' -R -o '+ wd+'/soap_out/'+key+'_soap"  ; } 2> '+key+'_soap.time &\n')
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
		f.write(' -R -o '+ wd+'/soap_out/'+key+'_soap"  ; } 2> '+key+'_soap.time &\n')
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

def prep_reads(n, length, pair): #filter length fq file and contig, and then combine
	f = open(wd+'/prep_reads.sh', 'w')
	job=0
	#merge .fq
	for (key, vdict) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		f.write(serverTag+' cat '+wd+'/fastq/'+key+'_1_sequence.txt ')
		if pair: f.write(wd+'/fastq/'+key+'_2_sequence.txt ')
		f.write('> '+wd+'/fastq/'+key+'.fq &\n')
		job+=1
		if job%nservers==0: f.write('wait\n')
	f.write('wait\n')
	# for (key, vdict) in seeds.items(): # copy sequence.txt to fq for Ray input requirement
		# serverTag = 'ssh '+servers[job%nservers]
		# f.write(serverTag+' cat '+wd+'/fastq/'+key+'_1_sequence.txt > '+wd+'/fastq/'+key+'_1.fq &\n')
		# if pair: f.write(serverTag+' cat '+wd+'/fastq/'+key+'_2_sequence.txt > '+wd+'/fastq/'+key+'_2.fq & \n')
		# job+=1
		# if job%nservers==0: f.write('wait\n')
	# f.write('wait\n')
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
		serverTag = 'ssh '+servers[job%nservers]
		i=1
		for fqfile in fqfiles:
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

def abyss(): #filter length fq file and contig, and then combine
	f = open(wd+'/abyss.sh', 'w')
	rval=[]
	job=0
	for (key, vdict) in seeds.items():
		serverTag = '{ time ssh '+servers[job%nservers]
		f.write(serverTag+' "cd '+wd+' && abyss-pe name='+key+' k='+abysskmer+' se='+wd+'/fastq/'+key+'abyss.fq" ; } 2> '+key+'_abyss.time &\n')
		rval.append(wd+'/'+key+'-unitigs.fa')
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
		f.write(dirscr+'partition.py '+wd+'/fastq/'+key+'.fq  100000\n')
	f.close()

def abyss_partition(): #filter length fq file and contig, and then combine
	f = open(wd+'/abyss_partition.sh', 'w')
	f1 = open(wd+'/abyss_combine.sh', 'w')
	rval=[]
	for (key, vdict) in seeds.items():
		rv=[]
		for j in xrange(20): #assume 100 chunks, unknown file size but some may not exist
			f.write('abyss-pe name='+key+'_chunk'+str(j)+' k='+abysskmer+' se='+wd+'/fastq/'+key+'.fq_'+str(j)+'  \n')
			rv.append(wd+'/'+key+'_chunk'+str(j)+'-unitigs.fa')
		f1.write('cat '+' '.join(rv) + ' > '+wd+'/'+key+'_abyss_partition.fa\n')
		rval.append(wd+'/'+key+'_abyss_partition.fa')
	f.close()
	return rval

def soap_partition(): #single end soap partition
	f = open(wd+'/soap_partition.sh', 'w')
	f1 = open(wd+'/soap_combine.sh', 'w')
	rval=[]
	for (key, vdict) in seeds.items():
		rv=[]
		for j in xrange(20): #assume 100 chunks, unknown file size but some may not exist
			f.write('abyss-pe name='+key+'_chunk'+str(j)+' k='+abysskmer+' se='+wd+'/fastq/'+key+'.fq_'+str(j)+'  \n')
			rv.append(wd+'/'+key+'_chunk'+str(j)+'-unitigs.fa')
		f1.write('cat '+' '.join(rv) + ' > '+wd+'/'+key+'_abyss_partition.fa\n')
		rval.append(wd+'/'+key+'_abyss_partition.fa')
	f.close()
	return rval

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
		f.write(' -R -o '+ wd+'/soap_out/'+key+'_soap"  ; } 2> '+key+'_soap.time &\n')
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
	
# { time ls; } 2> hello
def mira4(mira_local):
	f = open('mira.sh', 'w')
	job=0
	rval=[]
	for (key, vdict) in seeds.items():
		serverTag = '{ time ssh '+servers[job%nservers]
		ff = open(wd+'/'+key+'.conf', 'w')
		f.write(serverTag+' "cd '+wd+' && mira '+wd+'/'+key+'.conf" ; } 2> '+key+'_mira.time &\n')
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

def combineContig(contigLength1,r1,r2,r3,r4, r5, doMira, doPartition, pair):
# echo program key top1 top2 top3 N50 ncontigs n500 n1000 n2000 n5000 > contig.log
#abyssOut= wd+'/'+key+'-unitigs.fa'
#velvetOut=wd+'/velvet_'+key+'/contigs.fa'
# /mnt/cluster/xdeng/script/statContigs.py /mnt/san2/xdeng/Vhunt/131210_TN_CDC5/soap_out/Phan491DNA_soap.contig soap Phan491DNA
# /mnt/cluster/xdeng/script/statContigs.py /mnt/san2/xdeng/Vhunt/131210_TN_CDC5/velvet_Phan491DNA/contigs.fa meta_velvet Phan491DNA
# /mnt/cluster/xdeng/script/statContigs.py /mnt/san2/xdeng/Vhunt/131210_TN_CDC5/abyss_Phan491DNA/abyss-unitigs.fa abyss Phan491DNA
# /mnt/cluster/xdeng/script/statContigs.py Phan491DNA_assembly/Phan491DNA_d_results/Phan491DNA_out.unpadded.fasta mira Phan491DNA
# /mnt/cluster/xdeng/script/faLenFilter.py /mnt/san2/xdeng/Vhunt/131210_TN_CDC5/fastq/Phan491DNA_contig /mnt/san2/xdeng/Vhunt/131210_TN_CDC5/fastq/Phan491DNA_contig2 150
	f = open(wd+'/combineContig.sh', 'w')
	f1 = open(wd+'/statContig.sh', 'w')
	f2 = open(wd+'/blastContig.sh', 'w')
	f1.write('echo program key top1 top2 top3 N50 ncontigs n500 n1000 n2000 n5000\n')
	i=0
	for (key, vdict) in seeds.items():
		contig = wd+'/soap_out/'+key+'_soap.contig'
		if pair: 
			combfile=r1[i]+' '+r2[i]+' '+r3[i]
			if doMira:combfile+=' '+r4[i]
			if doPartition: combfile+=' '+r5[i]
		else:
			combfile=r1[i]
		f.write('cat '+combfile+ ' > '+wd+'/fastq/'+key+'_contig\n')
		f.write(dirscr+'faLenFilter.py '+wd+'/fastq/'+key+'_contig '+wd+'/fastq/'+key+'_contig2 '+str(contigLength1) +'\n')
		f1.write(dirscr+'statContigs.py '+r1[i]+'  soap '+key+'\n')
		f1.write(dirscr+'statContigs.py '+r2[i]+'  meta_velvet '+key+'\n')
		f1.write(dirscr+'statContigs.py '+r3[i]+'  abyss '+key+'\n')
		if doMira: f1.write(dirscr+'statContigs.py '+r4[i]+'  mira '+key+'\n')
		if doPartition: f1.write(dirscr+'statContigs.py '+r5[i]+'  partition '+key+'\n')
		f1.write(dirscr+'statContigs.py '+wd+'/fastq/'+key+'_contig3 ensemble '+key+'\n')
		f2.write(dirscr+'blastContig.py '+r1[i]+'  soap '+key+'\n')
		f2.write(dirscr+'blastContig.py '+r2[i]+'  meta_velvet '+key+'\n')
		f2.write(dirscr+'blastContig.py '+r3[i]+'  abyss '+key+'\n')
		if doMira: f2.write(dirscr+'blastContig.py '+r4[i]+'  mira '+key+'\n')
		if doPartition: f2.write(dirscr+'blastContig.py '+r5[i]+'  mira '+key+'\n')
		f2.write(dirscr+'blastContig.py '+wd+'/fastq/'+key+'_contig3 '+'  ensemble '+key+'\n')
		# /mnt/cluster/xdeng/script/blastContig.py /mnt/san2/xdeng/Vhunt/131210_TN_CDC5/soap_out/Phan491DNA_soap.contig soap Phan491DNA
# /mnt/cluster/xdeng/script/blastContig.py /mnt/san2/xdeng/Vhunt/131210_TN_CDC5/velvet_Phan491DNA/contigs.fa meta_velvet Phan491DNA
		i+=1
	f.close()
	f1.close()
	f2.close()

# def compareContigs(n, len, pair): #filter length fq file and contig, and then combine
	# f = open(wd+'/compareContigs.sh', 'w')
	# for (key, vdict) in seeds.items():
		# contig = wd+'/soap_out/'+key+'_soap.contig'
		# contig2 = wd+'/velvet_'+key+'/meta-velvetg.contigs.fa'
		# f.write(dirscr+'compareContigs.py '+contig+' '+contig2+'\n')
	# f.close()

def cap3():
#ssh bsidna3 cap3 /mnt/san2/xdeng/Vhunt/131210_TN_CDC5/fastq/Phan491DNA_contig2 &
	f = open('cap3.sh', 'w')
	job=0
	for (key, vdict) in seeds.items():
		serverTag = '{ time ssh '+servers[job%nservers]
		f.write(serverTag+' "cap3 '+wd+'/fastq/'+key+'_contig2" ; } 2> '+key+'_cap3.time &\n')
		job+=1
		if job%nservers==0: f.write('wait\n')
	f.close()

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
def combineContig_reads(n, skipread): #filter length fq file and contig, and then combine
	f = open(wd+'/combineContig_reads.sh', 'w')
	for (key, vdict) in seeds.items():
		contig = wd+'/fastq/'+key+'_contig4 ' #contig combined with reads
		combine = wd+'/fastq/'+key+'_c ' #contig combined with reads
		if skipread: f.write('cat '+contig+' > '+combine +'\n')
		else: f.write('cat '+contig+' '+wd+'/fastq/'+key+'.fa > '+combine +'\n')
		f.write(dirscr+'splitQuery.py '+combine+' '+str(n)+'\n')
	f.close()

def blastVirus(n, hsp):
	f = open(wd+'/blast_virus.sh', 'w')
	f1 = open(wd+'/blast_nr.sh', 'w')
	f2 = open(wd+'/blast_virus_parser.sh', 'w')
	f3 = open(wd+'/blast_nr_filter.sh', 'w')
	f4 = open(wd+'/blast_output_sort.sh', 'w')
	f5 = open(wd+'/blast_nr_mystery.sh', 'w')
	f6 = open(wd+'/blast_nr_mystery_filter.sh', 'w')
	f7 = open(wd+'/movetowww.sh', 'w')
	try: os.mkdir(wd+'/blast_virus_out/')
	except: pass
	try: os.mkdir(wd+'/blast_nr_out/')
	except: pass
	try: os.mkdir(wd+'/blast_filter_out/')
	except: pass
	try: os.mkdir(wd+'/'+base+'/')
	except: pass
	try: os.mkdir(wd+'/'+base+'/mystery/')
	except: pass
	try: os.mkdir(wd+'/'+base+'/hist/')
	except: pass
	try: os.mkdir(wd+'/blast_filter_out/tmp/')
	except: pass
	# f5 = open(wd+'/blast_filter_out/index.html', 'w')
	# f5.write('<html>')
	dbname= os.path.abspath(virusdbpath)
	dbnr= os.path.abspath(nvnrdbpath)
	job=0
	allout = 'cat '
	allmys = 'cat '
	allcombine='cat '
	for (key, vdict) in seeds.items():
		queryname = wd+'/fastq/'+key+'_c'
		mysfile=wd+'/'+base+'/mystery/'+key+'_m.fasta'
		outname = wd+'/blast_virus_out/'+key+'_blast.xml'
		outname1 = wd+'/blast_nr_out/'+key+'_blast.xml'
		outtxt = wd+'/blast_filter_out/'+key+'_blast_filter.txt'
		# htmlFile =os.path.basename(os.path.splitext(outtxt)[0]+'.html')
		# f5.write('<a href="'+htmlFile+'">'+htmlFile+'</a><br>')
		merge='cat '
		mys='cat '
		for i in xrange(n):
			serverTag = 'ssh '+servers[job%nservers]
			f.write(serverTag+' '+blastxpath+' -num_threads '+thread+'  -evalue 0.01 -outfmt 5 -dbsize 10000000 -searchsp 1000000000 -max_target_seqs 1')
			f.write(' -db_soft_mask 21 -best_hit_overhang 0.1 -best_hit_score_edge 0.1') #masking and best hit
			f1.write(serverTag+' '+blastxpath+' -num_threads '+thread+'  -evalue 0.01 -outfmt 5 -dbsize 10000000 -searchsp 1000000000 -max_target_seqs 1')
			f1.write(' -best_hit_overhang 0.1 -best_hit_score_edge 0.1') #best hit only
			subquery=queryname+'_'+str(i)
			sigfaname = subquery+'_s' #output significant virus fa file for NR
			submyspre = wd+'/blast_filter_out/tmp/'+key+'_prem.fasta'+'_'+str(i) #output mysterious contigs
			submys = wd+'/blast_filter_out/tmp/'+key+'_m.fasta'+'_'+str(i) #output mysterious contigs filtered by NR
			submysxml = wd+'/blast_filter_out/tmp/'+key+'_m.xml_'+str(i) #output mysterious contigs
			subxml=outname+'_'+str(i)
			subxml1=outname1+'_'+str(i)
			subout=outtxt+'_'+str(i)
			f.write(' -db '+dbname+ ' -query '+ subquery+' -out '+subxml+' & \n')
			f1.write(' -db '+dbnr+ ' -query '+ sigfaname+' -out '+subxml1+' & \n')
			f2.write(serverTag+' '+dirscr+'blast_parser.py '+subxml+ ' '+ subquery+ ' '+ sigfaname+ ' '+submyspre+' '+str(myslen)+' & \n')
			f5.write(serverTag+' '+blastxpath+' -num_threads '+thread+'  -evalue 0.001 -outfmt 5 -dbsize 10000000 -searchsp 1000000000 -max_target_seqs 1')
			f5.write(' -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -db '+dbnr+ ' -query '+ submyspre+' -out '+submysxml+' & \n')
			f6.write(serverTag+' '+dirscr+'blast_filter_NR_mys.py '+ submysxml +' '+submyspre+' '+ submys+' & \n')
			f3.write(serverTag+' '+dirscr+'blast_filter_NR.py '+ subxml + '  '+ subxml1+'  '+ subquery+' '+subout+' '+hsp+' & \n')
			job+=1
			if job%nservers==0: 
				f.write('wait\n')
				f1.write('wait\n')
				f2.write('wait\n')
				f3.write('wait\n')
				f5.write('wait\n')
				f6.write('wait\n')
			merge+=subout+' '
			mys+=submys+' '
			#rem+='rm '+submys+'\n'
		merge+='> '+outtxt
		mys+='> '+mysfile
		allout+= outtxt + ' '
		allmys+= mysfile + ' '
		f4.write(merge+'\n')
		f4.write(mys+'\n')
		#f4.write(rem+'\n')
		combine = wd+'/fastq/'+key+'.fa' 
		allcombine+= combine +' '
		f4.write(dirscr+'blast_output_sort.py '+ outtxt +' '+combine+'\n')
	outtxt = wd+'/blast_filter_out/'+'all'+'_blast_filter.txt'
	mysfile =wd+'/'+base+'/mystery/'+'all'+'_m.fasta'
	combine =wd+'/fastq/'+'all'+'.fa'
	allout +='> '+ outtxt
	allmys +='> '+ mysfile
	allcombine += '> '+ combine
	f4.write(allout+'\n')
	f4.write(allmys+'\n')
	f4.write(allcombine+'\n')
	f4.write(dirscr+'blast_output_sort.py '+ outtxt +' '+combine+'\n')
	f7.write('rm -rf '+wd+'/'+base+'/aln/\n')
	f7.write('rm -rf '+wd+'/'+base+'/fasta/\n')
	f7.write('rm -rf '+wd+'/'+base+'/pie/\n')
	f7.write('mv -f '+wd+'/blast_filter_out/aln '+wd+'/'+base+'/\n')
	f7.write('mv -f '+wd+'/blast_filter_out/fasta '+wd+'/'+base+'/\n')
	f7.write('mv -f '+wd+'/blast_filter_out/pie '+wd+'/'+base+'/\n')
	f7.write('mv -f '+wd+'/blast_filter_out/*.html '+wd+'/'+base+'/\n')
	f7.write('cp '+dirscr+'*.php '+wd+'/'+base+'/\n')
	f7.write('scp -r '+wd+'/'+base+' xdeng@dnasrv01dmzdr:/bsri/aspera_transfers/hold_test/\n')
	f.close()
	f1.close()
	f2.close()
	f3.close()
	f4.close()
	f5.close()
	f6.close()
	f7.close()
	
if __name__ == "__main__":
	fasta=False
	bam=False
	mira_local=False
	sra=False
	keep_human=False
	#keep_bac=False
	pair=True
	skipadaptor=False
	skipread=False
	doMira=False
	doPartition=False
	thread = '8'
	abysskmer='31'
	soapkmer='31'
	metakmer='31'
	hsp='NO' # show whole reads, if 'YES' show only blast hsp
	wd = os.path.abspath(os.path.dirname('.'))
	#genseedfile()
	# if pair:
	seeds=readSeeds2()
	# else:
		# seeds=readSeeds1()
	if not wd.startswith('/mnt/'): wd='/mnt'+wd
	base = os.path.basename(wd)
	print 'path',wd
	print 'thread', thread
	print 'base', base
	for key in seeds.keys():
		print key, seeds[key]
	print 'fasta', fasta
	print 'keep_human', keep_human
	# print 'keep_bac', keep_bac
	print 'pairend', pair
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
	# print >>oof, 'keep_bac', keep_bac
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
	length=100 #lengh reads for blast
	contigLength1 = 150 #length before cap3
	contigLength2 = 300 #length after cap3 before blast
	myslen = 1200
	print 'read threshold', length
	print 'contig before cap3', contigLength1
	print 'contig after cap3', contigLength2
	print 'mys threshold', myslen
	prep_sra()
	polyA()
	if fasta: fa2fq()
	if bam: bam2fq()
	# bowtieHuman(pair)
	# #bowtieNT()
	# skipbowtie()
	bowtieBac(pair, keep_human)
	sam2fq(pair, keep_human)

	check_pair_overlap()
	prepBlastFile_adaptor()
	trim()
	skip_adaptor()
	fq_check()
	prep_reads(n, length, pair)
	if pair: r1=soap_pair()
	else: r1=soap_single() 
	r2=meta_velvet()
	r3=abyss()
	r5=abyss_partition()
	r4=mira4(mira_local)
	# r6=Ray()
	combineContig(contigLength1,r1,r2,r3,r4,r5,doMira, doPartition, pair)
	cap3()
	cap3MergeSinglet(contigLength2)
	combineContig_reads(n, skipread)
	blastVirus(n, hsp)
	clean_dir()
	prepBlastFile()
	
	#if pair: prepPriceFile() #not doing price

	sf=open('pipeline_run.sh', 'w')
	if fasta: sf.write('source fa2fq.sh\nwait\n')
	if bam: sf.write('source bam2fq.sh\nwait\n')
	# os.system('PriceTI  -fp s1.fq s2.fq 300 -icf seeds.fa 1 1 5  -nc 30 -dbmax 72 -mol 35 -tol 20 -mpi 80 -target 90 2 2 2 -a 7 -o price.fa')

	if sra: 
		sf.write('source sra.sh\n')
	# if keep_human:
		# sf.write('source skipbowtie.sh\n')
	# else:
		# sf.write('source bowtieHuman.sh\nwait\n')
		# sf.write('source bowtiesam2fq.sh >>stats.log\nwait\n')

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
	#sf.write('source check_pair.sh >>stats.log \nwait\n')
	sf.write('source polyA_raw.sh >>stats.log \n')
	sf.write('source polyA_clean.sh >>stats.log \n')
	sf.write('source prep_reads.sh\nwait\n')
	sf.write('source soap.sh\nwait\n')
	#sf.write('source Ray.sh\n')
	sf.write('source meta_velvet.sh\nwait\n')
	sf.write('source abyss.sh\nwait\n')
	if doMira:
		sf.write('source mira.sh\nwait\n')
	if doPartition: 
		sf.write('source abyss_partition.sh\nwait\n')
		sf.write('source abyss_combine.sh\nwait\n')
	sf.write('source combineContig.sh \nwait\n')
	sf.write('source cap3.sh\nwait\n')
	sf.write('source cap3MergeSinglet.sh \nwait\n')
	sf.write('source statContig.sh > contig.log\n')
	sf.write('source blastContig.sh > contig.blast\n')
	sf.write('source combineContig_reads.sh >>stats.log \nwait\n')
	sf.write('source blast_virus.sh\nwait\n')
	sf.write('source blast_virus_parser.sh >>stats.log  \nwait\n')
	sf.write('source blast_nr_mystery.sh \nwait\n')
	sf.write('source blast_nr_mystery_filter.sh >>stats.log  \nwait\n')
	sf.write('source blast_nr.sh\nwait\n')
	sf.write('source blast_nr_filter.sh >>stats.log  \nwait\n')
	sf.write('source blast_output_sort.sh >>stats.log \nwait\n')
	sf.write('firstpage.py\n')
	sf.write('source prepBlastFile.sh \n')
	# if pair: sf.write('source prepPriceFile.sh \n')
	sf.write('source movetowww.sh \n')
	sf.close()
