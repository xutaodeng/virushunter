#!/usr/bin/env python
from collections import defaultdict
import os.path
import os
import string
import sys

soappath='/home/ec2-user/mnt/tools/SOAPdenovo2-src-r239/SOAPdenovo-63mer'
bowtiepath='/home/ec2-user/mnt/tools/bowtie2-2.1.0/bowtie2' 
blastxpath='/home/ec2-user/mnt/tools/ncbi-blast-2.2.28+/bin/blastx'
bowtieindexpath='/home/ec2-user/mnt/bowtieIndex/mrnadna_bowtie'
virusdbpath='/home/ec2-user/mnt/virus/virus_mask'
nvnrdbpath='/home/ec2-user/mnt/refseq/nvrefseq'
dirscr='/home/ec2-user/mnt/script/'

#$parameter = $fq1.' '.$fq2.' '.$adaptor.' '.$trim.' '.$qtrim.' '.$platform.' '.$pair.' '.$human.' '.$contig;

fq1 = sys.argv[1]
fq2 = sys.argv[2]
adaptors = sys.argv[3] #dirscr+'adaptor.fa'
skipadaptor=sys.argv[4] #=False
qtrim = sys.argv[5]
platform=sys.argv[6] #=True
pair=sys.argv[7]
skiphuman=sys.argv[8] #True
contig=sys.argv[9] # read, contig, both

if platform=="ilmn13": platform= '64'
elif platform=="ilmn15": platform= '64'
elif platform=="ilmn18": platform= '33'
elif platform=="notsure": platform= '33'

if qtrim=="phred7": qtrim= '7'
elif qtrim=="phred10": qtrim= '10'
elif qtrim=="phred13": qtrim= '13'
elif qtrim=="phred20": qtrim= '20'
elif qtrim=="no": qtrim= '0'

if skiphuman =='noremove': skiphuman=False
else: skiphuman=True
if pair=='pair-end': pair=True
else: pair=False
if skipadaptor =='trim': skipadaptor=True
else: skipadaptor=False

seeds=defaultdict(list)
key='library'
seeds[key].append(fq1)
seeds[key].append(fq2)
print 'skiphuman', skiphuman
print 'pairend', pair
print 'skipadaptor', skipadaptor

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
	
def skip_adaptor():
	f = open('skipadaptor.sh', 'w')
	for (key, fqfiles) in seeds.items():
		i=1
		for fqfile in fqfiles:
		#ff = open(wd+'/soap_config/'+key+'_soap.config', 'w')
			index=key+'_'+str(i)
			fid = wd+'/fastq/'+index
			f.write('recodeID.py '+fid+'.dup '+' '+wd+'/fastq/' + key+'_'+str(i)+'_sequence.txt '+key+ ' '+str(i)+'\n') #quality trimming
	f.close()

def prepBlastFile_adaptor(): #for adaptor filtering
	try: os.mkdir(wd+'/blast_filter_out/')
	except: pass
	try: os.mkdir(wd+'/blast_filter_out/blast2/')
	except: pass
	f = open('prepBlastFile_adaptor.sh', 'w')
	path='fastq/'
	f.write('cd fastq/ \n')
	bpath=wd+'/blast_filter_out/blast2/'
	directory = os.path.basename(wd)
	alldb=[]
	for (key, fqfiles) in seeds.items():
		i=1
		pairdb=[]
		for fqfile in fqfiles:
			index=key+'_'+str(i)
			fqfil=index+'.fil'
			f.write('makeblastdb -dbtype nucl -parse_seqids -in '+index+'.fa -out '+index+'_blastdb\n')
			pairdb.append(index+'_blastdb')
			i+=1
		f.write('blastdb_aliastool -dblist '+'\"'+' '.join(pairdb)+'\" -dbtype nucl -out '+key+'_blastdb -title \"'+key+'_blastdb\"\n')
		alldb.append(key+'_blastdb')
	f.write('blastdb_aliastool -dblist '+'\"'+' '.join(alldb)+'\" -dbtype nucl -out '+directory+'_blastdb -title \"'+directory+'_blastdb\"\n')
	f.write('mv *_blastdb* '+bpath+'\n')
	f.write('cd ..\n')
	f.close()

def trim():
	f = open('clonetrim.sh', 'w')
	f1 = open('clonefq2fa.sh', 'w')
	f2 = open('blast_adaptor.sh', 'w')
	f3 = open('blasttrim.sh', 'w')
	f4 = open('qualitytrim.sh', 'w')
	bpath=wd+'/blast_filter_out/blast2/'
	job=0
	for (key, fqfiles) in seeds.items():
		i=1
		for fqfile in fqfiles:
		#ff = open(wd+'/soap_config/'+key+'_soap.config', 'w')
			index=key+'_'+str(i)
			fid = wd+'/fastq/'+index
			fastqfile=wd+'/fastq/'+fqfile
			fqfil=wd+'/fastq/'+index+'.fil'
			#serverTag = 'ssh '+servers[job%nservers]
			serverTag = ''
			job+=1
			f.write(serverTag + ' '+dirscr+'clone_reads_rm_pair.py '+fqfil+' '+fid+'.dup \n')
			f1.write(dirscr+'fq2faID.py '+fid+'.dup '+' '+index+' '+fid+'.fa\n')
			f2.write(serverTag + ' blastn -task blastn -evalue 1  -max_target_seqs 100000000 -outfmt \'"6  qseqid  sseqid evalue qstart qend sstart send"\'')
			f2.write(' -query '+adaptors+ ' -num_threads 7 -db '+bpath+index+'_blastdb'+' -out '+fid+'.tab  \n')
			f3.write(serverTag + ' '+dirscr+'blast_trim.py '+fid+'.dup '+fid+'.tab '+fid+'.fq '+ key+ ' '+str(i)+' \n') #adaptor trimming
			f4.write(serverTag + ' '+dirscr+'trim_quality.py '+fid+'.fq '+wd+'/fastq/' + key+'_'+str(i)+'_sequence.txt '+platform+' '+fastqfile+' '+qtrim+' \n') #quality trimming
			i+=1
	f.close()
	f1.close()
	f2.close()
	f3.close()
	f4.close()

def check_pair_overlap():
	f = open('check_pair.sh', 'w')
	job=0
	for (key, fqfiles) in seeds.items():
		#serverTag = 'ssh '+servers[job%nservers]
		serverTag = ''
		job+=1
		f.write(serverTag + ' '+dirscr+'check_pair.py '+wd+'/fastq/' + key+'_1_sequence.txt '+wd+'/fastq/' + key+'_2_sequence.txt '+' \n') #quality trimming
	f.close()


def clean_dir():
	f2 = open('clean.sh', 'w')
	f2.write('find . -type f ! -name "*.gz"  -a ! -name "*.fastq" -delete\n')
	f2.close()

def polyA():
	f = open('polyA_raw.sh', 'w')
	f2 = open('polyA_clean.sh', 'w')
	for (key, fqfiles) in seeds.items():
		i=1
		for fqfile in fqfiles:
			f.write(dirscr+'polyA.py fastq/'+fqfile+' raw\n')
			f2.write(dirscr+'polyA.py fastq/' + key+'_'+str(i)+'_sequence.txt\n')
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
		i=1
		pairfa=[]
		for fqfile in fqfiles:
			index=key+'_'+str(i)
			f.write(dirscr+'fq2fa.py '+index +'_sequence.txt '+index+' '+index+'.fa 20 \n')
			pairfa.append(index+'.fa')
			i+=1
		f.write('cat '+' '.join(pairfa)+' > '+key+'.fa\n')
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
				# f.write(serverTag+' '+bowtiepath+' --no-hd --quiet --reorder --local -p 8 -x /mnt/cluster/xdeng/nt/nt'+str(j)+ \
				# ' -U '+wd+'/fastq/'+index +'_sequence.txt.tmp -S '+wd+'/fastq/'+ index+'_'+str(j)+'.sam \n')
				# job+=1
				# if job%nservers==0: f.write('wait\n')
			# f2.write(dirscr+'sam2fq.py '+wd+'/fastq/' + index+' '+wd+'/fastq/'+index+'_sequence.txt \n')
			# i+=1
	# f.close()
	# f2.close()

def bowtieHuman():
	f = open('bowtieHuman.sh', 'w')
	f2 = open('bowtiesam2fq.sh', 'w')
	job=0
	job2=0
	for (key, fqfiles) in seeds.items():
		i=1
		sams=[]
		fqfils=[]
		for fqfile in fqfiles:
			fastqfile=wd+'/fastq/'+fqfile
			index=key+'_'+str(i)
			sam=wd+'/fastq/'+ index+'.sam'
			sams.append(sam)
			fqfil=wd+'/fastq/'+index+'.fil'
			fqfils.append(fqfil)
			#serverTag = 'ssh '+servers[job%nservers]
			serverTag = ''
			f.write(serverTag+' '+bowtiepath+' --quiet --local --no-hd --reorder -p 7 -x '+bowtieindexpath+' -U '+fastqfile+' -S '+sam+' \n')
			job+=1
			#if job%nservers==0: f.write('wait\n')
			# f.write('/mnt/cluster/tools/bwa-0.5.9/bwa aln -n 0.2 -o 3 -t 7 /mnt/cluster/xdeng/hg19/mrnadna.fa fastq/' +index +'_sequence.txt.tmp > fastq/'+ index+'.sai \n')
			# f.write('/mnt/cluster/tools/bwa-0.5.9/bwa samse /mnt/cluster/xdeng/hg19/mrnadna.fa fastq/' + index+'.sai fastq/'+ index+'_sequence.txt.tmp > fastq/'+ index+'.sam \n')
			i+=1
		#serverTag2 = 'ssh '+servers[job2%nservers]
		serverTag2 = ''
		if pair: f2.write(serverTag2 + ' '+dirscr+'sam2fq.py '+sams[0]+' '+sams[1]+' '+fqfils[0]+' '+fqfils[1]+ '  \n')
		else: f2.write(serverTag2 + ' '+dirscr+'sam2fq.py '+sams[0]+' '+fqfils[0]+ '  \n')
		job2+=1
		#if job2%nservers==0: f2.write('wait\n')
	f.close()
	f2.close()
	

# now generate commandline sh for SOAPDENOVO
def soap_single():
	f = open('soap_single.sh', 'w')
	try: os.mkdir(wd+'/soap_config/')
	except: pass
	try: os.mkdir(wd+'/soap_out/')
	except: pass
	job=0
	for (key, vdict) in seeds.items():
		#serverTag = 'ssh '+servers[job%nservers]
		serverTag = ''
		ff = open(wd+'/soap_config/'+key+'_soap.config', 'w')
		f.write(serverTag+' '+soappath+' all -K 63')
		f.write(' -s '+ wd+'/soap_config/'+key+'_soap.config ')
		f.write(' -R -o '+ wd+'/soap_out/'+key+'_soap  \n')
		ff.write('max_rd_len=300\n[LIB]\n')
		ff.write('nreverse_seq=0\nasm_flags=3\nrank=1\n')
		ff.write('q='+wd+'/fastq/'+key+'_1_sequence.txt\n')
		ff.close()
		job+=1
		#if job%nservers==0: f.write('wait\n')
	f.close()
	
def soap_pair():
	f = open('soap_pair.sh', 'w')
	try: os.mkdir(wd+'/soap_config/')
	except: pass
	try: os.mkdir(wd+'/soap_out/')
	except: pass
	job=0
	for (key, vdict) in seeds.items():
		#serverTag = 'ssh '+servers[job%nservers]
		serverTag = ''
		ff = open(wd+'/soap_config/'+key+'_soap.config', 'w')
		f.write(serverTag+' '+soappath+' all -K 63')
		f.write(' -s '+ wd+'/soap_config/'+key+'_soap.config ')
		f.write(' -R -o '+ wd+'/soap_out/'+key+'_soap  \n')
		ff.write('max_rd_len=300\n[LIB]\n')
		ff.write('avg_ins=350\nreverse_seq=0\nasm_flags=3\nrank=1\n')
		ff.write('q1='+wd+'/fastq/'+key+'_1_sequence.txt\n'+'q2='+wd+'/fastq/'+key+'_2_sequence.txt \n')
		ff.close()
		job+=1
		#if job%nservers==0: f.write('wait\n')
	f.close()

def skip_read(n, len): #filter length fq file and contig, and then combine
	f = open(wd+'/skip_read.sh', 'w')
	for (key, vdict) in seeds.items():
		contig = wd+'/soap_out/'+key+'_soap.contig'
		contigfilter = wd+'/fastq/'+key+'_f'
		combine = wd+'/fastq/'+key+'_c' #contig combined with reads
		f.write(dirscr+'faLenFilter.py '+contig+' '+contigfilter+' '+str(len)+'\n')
		f.write('cp '+contigfilter+' '+combine +'\n')
		f.write(dirscr+'splitQuery.py '+combine+' '+str(n)+'\n')
	f.close()

def skip_contig(n, len): #filter length fq file and contig, and then combine
	f = open(wd+'/skip_contig.sh', 'w')
	for (key, vdict) in seeds.items():
		f.write('cat '+wd+'/fastq/'+key+'_1_sequence.txt ')
		if pair: f.write(wd+'/fastq/'+key+'_2_sequence.txt ')
		f.write(' > '+wd+'/fastq/'+key+'.fq \n')
		f.write(dirscr+'fq2fa.py ')
		f.write(wd+'/fastq/'+key+'.fq '+wd+'/fastq/'+key+'.fa '+str(len)+' \n')
		combine = wd+'/fastq/'+key+'_c' #contig combined with reads
		f.write('cp '+wd+'/fastq/'+key+'.fa '+combine +'\n')
		f.write(dirscr+'splitQuery.py '+combine+' '+str(n)+'\n')
	f.close()

def combineContig_reads(n, len, pair): #filter length fq file and contig, and then combine
	f = open(wd+'/combineContig_reads.sh', 'w')
	for (key, vdict) in seeds.items():
		f.write('cat '+wd+'/fastq/'+key+'_1_sequence.txt ')
		if pair: f.write(wd+'/fastq/'+key+'_2_sequence.txt ')
		f.write(' > '+wd+'/fastq/'+key+'.fq \n')
		f.write(dirscr+'fq2fa.py ')
		f.write(wd+'/fastq/'+key+'.fq '+wd+'/fastq/'+key+'.fa '+str(len)+' \n')
		contig = wd+'/soap_out/'+key+'_soap.contig'
		contigfilter = wd+'/fastq/'+key+'_f'
		combine = wd+'/fastq/'+key+'_c' #contig combined with reads
		f.write(dirscr+'faLenFilter.py '+contig+' '+contigfilter+' '+str(len)+'\n')
		f.write('cat '+contigfilter+' '+wd+'/fastq/'+key+'.fa > '+combine +'\n')
		f.write(dirscr+'splitQuery.py '+combine+' '+str(n)+'\n')
	f.close()
	
def fq2fa_all_reads(pair): #get all reads for getting paired read _r
	f = open(wd+'/fq2fa_all_reads.sh', 'w')
	for (key, vdict) in seeds.items():
		f.write('cat '+wd+'/fastq/'+key+'_1_sequence.txt ')
		if pair: f.write(wd+'/fastq/'+key+'_2_sequence.txt ')
		f.write(' > '+wd+'/fastq/'+key+'.fq \n')
		f.write(dirscr+'fq2fa.py ')
		f.write(wd+'/fastq/'+key+'.fq '+wd+'/fastq/'+key+'_r 0 '+' \n')
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
			#serverTag = 'ssh '+servers[job%nservers]
			serverTag = ''
			f.write(serverTag+' '+blastxpath+' -num_threads 8 -evalue 0.01 -outfmt 5 -dbsize 10000000 -searchsp 1000000000 -max_target_seqs 1')
			f.write(' -db_soft_mask 21 -best_hit_overhang 0.1 -best_hit_score_edge 0.1') #masking and best hit
			f1.write(serverTag+' '+blastxpath+' -num_threads 8 -evalue 0.01 -outfmt 5 -dbsize 10000000 -searchsp 1000000000 -max_target_seqs 1')
			f1.write(' -best_hit_overhang 0.1 -best_hit_score_edge 0.1') #best hit only
			subquery=queryname+'_'+str(i)
			sigfaname = subquery+'_s' #output significant virus fa file for NR
			submyspre = wd+'/blast_filter_out/tmp/'+key+'_prem.fasta'+'_'+str(i) #output mysterious contigs
			submys = wd+'/blast_filter_out/tmp/'+key+'_m.fasta'+'_'+str(i) #output mysterious contigs filtered by NR
			submysxml = wd+'/blast_filter_out/tmp/'+key+'_m.xml_'+str(i) #output mysterious contigs
			subxml=outname+'_'+str(i)
			subxml1=outname1+'_'+str(i)
			subout=outtxt+'_'+str(i)
			f.write(' -db '+dbname+ ' -query '+ subquery+' -out '+subxml+'  \n')
			f1.write(' -db '+dbnr+ ' -query '+ sigfaname+' -out '+subxml1+'  \n')
			f2.write(serverTag+' '+dirscr+'blast_parser.py '+subxml+ ' '+ subquery+ ' '+ sigfaname+ ' '+submyspre+'  \n')
			f5.write(serverTag+' '+blastxpath+' -num_threads 8 -evalue 0.001 -outfmt 5 -dbsize 10000000 -searchsp 1000000000 -max_target_seqs 1')
			f5.write(' -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -db '+dbnr+ ' -query '+ submyspre+' -out '+submysxml+'  \n')
			f6.write(serverTag+' '+dirscr+'blast_filter_NR_mys.py '+ submysxml +' '+submyspre+' '+ submys+'  \n')
			f3.write(serverTag+' '+dirscr+'blast_filter_NR.py '+ subxml + '  '+ subxml1+'  '+ subquery+' '+subout+' '+hsp+'  \n')
			job+=1
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
		combine = wd+'/fastq/'+key+'_r' 
		allcombine+= combine +' '
		f4.write(dirscr+'blast_output_sort.py '+ outtxt +' '+combine+'\n')
	outtxt = wd+'/blast_filter_out/'+'all'+'_blast_filter.txt'
	mysfile =wd+'/'+base+'/mystery/'+'all'+'_m.fasta'
	combine =wd+'/fastq/'+'all'+'_r'
	allout +='> '+ outtxt
	allmys +='> '+ mysfile
	allcombine += '> '+ combine
	f4.write(allout+'\n')
	f4.write(allmys+'\n')
	f4.write(allcombine+'\n')
	f4.write(dirscr+'blast_output_sort.py '+ outtxt +' '+combine+'\n')
	f7.write('mv '+wd+'/blast_filter_out/aln '+wd+'/'+base+'/\n')
	f7.write('mv '+wd+'/blast_filter_out/fasta '+wd+'/'+base+'/\n')
	f7.write('mv '+wd+'/blast_filter_out/mystery '+wd+'/'+base+'/\n')
	f7.write('mv '+wd+'/blast_filter_out/pie '+wd+'/'+base+'/\n')
	f7.write('mv '+wd+'/blast_filter_out/*.html '+wd+'/'+base+'/\n')
	f7.write('cp '+dirscr+'*.php '+wd+'/'+base+'/\n')
	f.close()
	f1.close()
	f2.close()
	f3.close()
	f4.close()
	f5.close()
	f6.close()
	f7.close()

if __name__ == "__main__":
	hsp='NO' # show whole reads, if 'YES' show only blast hsp
	#print 'hsp', hsp
	wd = os.path.abspath(os.path.dirname('.'))
	if not wd.startswith('/mnt/'): wd='/mnt'+wd
	base = os.path.basename(wd)
	print 'path',wd
	print 'base', base
	for key in seeds.keys():
		print key, seeds[key]
	n=60
	len=60 #len of query (contig, and reads) threadhold

	polyA()
	bowtieHuman()
	#bowtieNT()
	skipbowtie()
	check_pair_overlap()
	prepBlastFile_adaptor()
	trim()
	skip_adaptor()
	if pair: soap_pair()
	else: soap_single() 
	skip_read(n, len)
	skip_contig(n, len)
	combineContig_reads(n, len, pair)
	blastVirus(n, hsp)
	fq2fa_all_reads(pair)
	clean_dir()
	prepBlastFile()
	if pair: prepPriceFile()
	sf=open('pipeline_run.sh', 'w')
	# os.system('PriceTI  -fp s1.fq s2.fq 300 -icf seeds.fa 1 1 5  -nc 30 -dbmax 72 -mol 35 -tol 20 -mpi 80 -target 90 2 2 2 -a 7 -o price.fa')

	if skiphuman:
		sf.write('echo "skipping human removal.." >> run.log \n')
		sf.write('source skipbowtie.sh\n')
	else:
		sf.write('echo "removing human reads.." >> run.log\n')
		sf.write('source bowtieHuman.sh\n')
		sf.write('source bowtiesam2fq.sh >>stats.log\n')
	sf.write('source clonetrim.sh >>stats.log \n')
	
	if skipadaptor:
		sf.write('echo "skipping trimming.." >> run.log\n')
		sf.write('source skipadaptor.sh\n')
	else:
		sf.write('echo "trimming reads.." >> run.log\n')
		sf.write('source clonefq2fa.sh \n')
		sf.write('source prepBlastFile_adaptor.sh \n')
		sf.write('source blast_adaptor.sh \n')
		sf.write('source blasttrim.sh >>stats.log \n')
		sf.write('source qualitytrim.sh >>stats.log \n')
	sf.write('echo "checking paired-end overlap.." >> run.log\n')
	sf.write('source check_pair.sh >>stats.log \n')
	sf.write('source polyA_raw.sh >>stats.log \n')
	sf.write('source polyA_clean.sh >>stats.log \n')
	if contig!='read' and pair:
		sf.write('echo "SOAP2 denovo assembly.." >> run.log\n')
		sf.write('source soap_pair.sh\n')
	elif contig!='read' and not pair:
		sf.write('echo "SOAP2 denovo assembly.." >> run.log\n')
		sf.write('source soap_single.sh\n')

	if contig=='contig':
		sf.write('source skip_read.sh >>stats.log \n')
	elif contig=='read': 
		sf.write('source skip_contig.sh >>stats.log \n')
	else:
		sf.write('source combineContig_reads.sh >>stats.log \n')
	sf.write('echo "BLASTX virus.." >> run.log\n')
	sf.write('source blast_virus.sh\n')
	sf.write('source blast_virus_parser.sh >>stats.log  \n')
	sf.write('source blast_nr_mystery.sh \n')
	sf.write('source blast_nr_mystery_filter.sh >>stats.log  \n')
	sf.write('echo "BLASTX non-vurus NR.." >> run.log\n')
	sf.write('source blast_nr.sh\n')
	sf.write('source blast_nr_filter.sh >>stats.log  \n')
	sf.write('echo "Finalizing.." >> run.log\n')
	sf.write('source fq2fa_all_reads.sh\n')
	sf.write('source blast_output_sort.sh >>stats.log \n')
	sf.write('firstpage.py\n')
	sf.write('source prepBlastFile.sh \n')
	if pair: sf.write('source prepPriceFile.sh \n')
	sf.write('source movetowww.sh \n')
	sf.close()
