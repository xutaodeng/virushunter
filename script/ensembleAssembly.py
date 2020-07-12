#!/usr/bin/env python
from collections import defaultdict
import os.path
import os
import string
import sys

soappath='/mnt/cluster/tools/SOAPdenovo-63mer'
velvetg='/mnt/cluster/tools/velvet_1.2.10/velvetg'
velveth='/mnt/cluster/tools/velvet_1.2.10/velveth'
meta_velvetg="/mnt/cluster/tools/MetaVelvet-1.2.02/meta-velvetg"
Minimo="/mnt/cluster/tools/amos-3.1.0/bin/Minimo"
abysspath='/mnt/cluster/tools/abyss/bin/abyss-pe'
cap3path='/mnt/cluster/tools/CAP3/cap3'
dirscr='/mnt/cluster/xdeng/script/'

# PE= pe 180 20 /FULL_PATH/frag_1.fastq /FULL_PATH/frag_2.fastq
# SE= FULL_PATH/short_1.fastq
# OUTDIR=FULL_PATH_OUTPUTDIR/
# NUM_THREADS= 16
# SOAP_KMER=31
# ABYSS_KMER = 31
# METAVELVET_KMER=31
# CON_LEN_DBG=150
# CON_LEN_OLC=300
# ASSEMBLY_MODE='quick', 'optimal'

def readConfig(configfile):
	global pair, thread, abysskmer, soapkmer, metakmer, insert_size, insert_sd, seeds, contigLength1, contigLength2, mode
	seeds=defaultdict(list)
	f=open(configfile, 'r')
	for line in f:
		if line.strip().startswith('#'): continue
		parts = line.strip().split('=',1)
		name, value = parts[0].strip(), parts[1].strip()
		if name == 'PE': 
			pair=True
			insert_size, insert_sd, file1, file2 = value.split()
			seeds['ensemble'].append(file1)
			seeds['ensemble'].append(file2)
		if name == 'SE': 
			pair=False
			file1 = value
			seeds['ensemble'].append(file1)
		if name == 'NUM_THREADS': thread = value
		if name == 'ABYSS_KMER': abysskmer=value
		if name == 'SOAP_KMER':soapkmer=value
		if name == 'METAVELVET_KMER':metakmer=value
		if name == 'CON_LEN_DBG': contigLength1=int(value)
		if name == 'CON_LEN_OLC': contigLength2=int(value)
		if name == 'ASSEMBLY_MODE': mode=value
	f.close()

def prep_reads(pair): #filter length fq file and contig, and then combine
	f = open(wd+'/prep_reads.sh', 'w')
	try: os.mkdir(wd+'/fastq/')
	except: pass
	#below prep abyss
	for (key, fastqs) in seeds.items():
		f.write('cat '+fastqs[0]+' ')
		if pair: f.write(fastqs[1])
		f.write(' > '+wd+'/fastq/'+key+'.fq\n')
	for (key, vdict) in seeds.items():
		f.write(dirscr+'fqLenFilter.py '+ wd+'/fastq/'+key+'.fq '+' '+wd+'/fastq/'+key+'abyss.fq 35\n')
	f.close()

def meta_velvet():
	f = open('meta_velvet.sh', 'w')
	# velveth out-dir 51 -fastq -shortPaired HMP.small/SRR041654_shuffled.fastq  HMP.small/SRR041655_shuffled.fastq 
	# velvetg out-dir -exp_cov auto -ins_length 260
	# meta-velvetg out-dir -ins_length 260 | tee logfile
	rval=[]
	for (key, fastqs) in seeds.items():
		outdir = wd+'/velvet_'+key+'/'
		if pair: 
			f.write(velveth+' '+outdir+' '+metakmer)
			f.write(' -shortPaired -fastq '+fastqs[0]+' '+fastqs[1]+' && ')
			f.write(velvetg+' '+outdir+' -exp_cov auto -ins_length '+insert_size+' && ')
			f.write(meta_velvetg+' '+outdir+' -ins_length '+insert_size+' \n')
		else:
			f.write(velveth+' '+outdir+' '+metakmer)
			f.write(' -short -fastq '+fastqs[0]+' && ')
			f.write(velvetg+' '+outdir+' -exp_cov auto && ')
			f.write(meta_velvetg+' '+outdir+'\n')
		velvetOut=wd+'/velvet_'+key+'/contigs.fa'
		rval.append(velvetOut)
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
	for (key, fastqs) in seeds.items():
		ff = open(wd+'/soap_config/'+key+'_soap.config', 'w')
		f.write(soappath+' all -K '+ soapkmer)
		f.write(' -s '+ wd+'/soap_config/'+key+'_soap.config ')
		f.write(' -R -o '+ wd+'/soap_out/'+key+'_soap\n')
		ff.write('max_rd_len=300\n[LIB]\n')
		ff.write('nreverse_seq=0\nasm_flags=3\nrank=1\n')
		ff.write('q='+fastqs[0]+'\n')
		ff.close()
		soapOut=wd+'/soap_out/'+key+'_soap.contig'
		rval.append(soapOut)
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
	for (key, fastqs) in seeds.items():
		ff = open(wd+'/soap_config/'+key+'_soap.config', 'w')
		f.write(soappath+' all -K '+soapkmer)
		f.write(' -s '+ wd+'/soap_config/'+key+'_soap.config ')
		f.write(' -R -o '+ wd+'/soap_out/'+key+'_soap\n')
		ff.write('max_rd_len=300\n[LIB]\n')
		ff.write('avg_ins='+insert_size+'\nreverse_seq=0\nasm_flags=3\nrank=1\n')
		ff.write('q1='+fastqs[0]+'\n'+'q2='+fastqs[1]+' \n')
		ff.close()
		soapOut=wd+'/soap_out/'+key+'_soap.contig'
		rval.append(soapOut)
	f.close()
	return rval

def abyss(): #filter length fq file and contig, and then combine
	f = open(wd+'/abyss.sh', 'w')
	rval=[]
	job=0
	for (key, vdict) in seeds.items():
		try: os.mkdir(wd+'/abyss_'+key)
		except: pass
		f.write('cd '+wd+' && abyss-pe -C '+wd+'/abyss_'+key+' name='+key+' k='+abysskmer+' se='+wd+'/fastq/'+key+'abyss.fq\n')
		rval.append(wd+'/abyss_'+key+'/'+key+'-unitigs.fa')
	f.close()
	return rval

def partition(): #filter length fq file and contig, and then combine
	f = open(wd+'/partition.sh', 'w')
	for (key, vdict) in seeds.items():
		f.write(dirscr+'partition.py '+wd+'/fastq/'+key+'abyss.fq  100000\n')
	f.close()

def abyss_partition(): #filter length fq file and contig, and then combine
	f = open(wd+'/abyss_partition.sh', 'w')
	f1 = open(wd+'/abyss_combine.sh', 'w')
	rval=[]
	for (key, vdict) in seeds.items():
		rv=[]
		for j in xrange(20): #assume 100 chunks, unknown file size but some may not exist
			try: os.mkdir(wd+'/abyss_'+key+'_chunk'+str(j))
			except: pass
			f.write('cd '+wd+' && abyss-pe -C '+wd+'/abyss_'+key+'_chunk'+str(j)+' name='+key+'_chunk'+str(j)+ \
			' k='+abysskmer+' se='+wd+'/fastq/'+key+'abyss.fq_'+str(j)+'\n')
			#f.write('abyss-pe name='+key+'_chunk'+str(j)+' k='+abysskmer+' se='+wd+'/fastq/'+key+'.fq_'+str(j)+'  \n')
			outfile = wd+'/abyss_'+key+'_chunk'+str(j)+'/'+key+'_chunk'+str(j)+'-unitigs.fa'
			rv.append(outfile)
		f1.write('cat '+' '.join(rv) + ' > '+wd+'/'+key+'_abyss_partition.fa\n')
		rval.append(wd+'/'+key+'_abyss_partition.fa')
	f.close()
	f1.close()
	return rval

def combineContig2(contigLength1,assembly_para):
	f = open(wd+'/combineContig2.sh', 'w')
	combfile=[]
	contig_out = wd+'/contig_individual/'
	for (key, vdict) in seeds.items():
		if 'S' in assembly_para: combfile.append(contig_out+key+'_S.contig') #soap
		if 'V' in assembly_para: combfile.append(contig_out+key+'_V.contig') #velvet
		if 'A' in assembly_para: combfile.append(contig_out+key+'_A.contig') 
		if 'a' in assembly_para: combfile.append(contig_out+key+'_a.contig') 
		f.write('cat '+' '.join(combfile)+ ' > '+wd+'/fastq/'+key+'_contig\n')
		combfile=[]
		f.write(dirscr+'faLenFilter.py '+wd+'/fastq/'+key+'_contig '+wd+'/fastq/'+key+'_contig2 '+str(contigLength1) +'\n')
	f.close()

#move contig for individual
def moveContig1(assembly_para, r1,r2,r3,r4, r5,r6,r7, r8, r9, r10, r11):
	contig_out = wd+'/contig_individual/'
	try: os.mkdir(contig_out)
	except: pass
	f = open(wd+'/moveContig.sh', 'w')
	i = 0
	for (key, vdict) in seeds.items():
		if 'S' in assembly_para: f.write('cp '+r1[i]+' '+contig_out+key+'_S.contig\n') #soap
		if 'V' in assembly_para:  f.write('cp '+r2[i]+' '+contig_out+key+'_V.contig\n') #velvet
		if 'A' in assembly_para: 
			f.write('cp --dereference '+r3[i]+' '+contig_out+key+'_A.contig\n')
		if 'a' in assembly_para: f.write('cp '+r7[i]+' '+contig_out+key+'_a.contig\n')
		i+=1
	f.close()

#move contig for ensembled assemblers
def moveContig2(assembly_para):
	contig_out = wd+'/contig_'+assembly_para+'/'
	try: os.mkdir(contig_out)
	except: pass
	f = open(contig_out+'moveContigIndividual.sh', 'w')
	i = 0
	for (key, vdict) in seeds.items():
		key2=key+'OLC'
		if 'C' in assembly_para: f.write('cp '+wd+'/fastq/'+key+'_contig3 '+contig_out+key+'_C.contig\n')
		if 'O' in assembly_para: f.write('cp '+wd+'/fastq/'+key+'_contig2-contigs.fa '+contig_out+key+'_O.contig\n')
		i+=1
	f.close()

def cap3():
	f = open('cap3.sh', 'w')
	for (key, vdict) in seeds.items():
		f.write(cap3path+' '+wd+'/fastq/'+key+'_contig2\n')
	f.close()

def minimo():
	f = open('minimo.sh', 'w')
	for (key, vdict) in seeds.items():
		f.write(Minimo+' '+wd+'/fastq/'+key+'_contig2 -D FASTA_EXP=1\n')
	f.close()

def minimoFilter(contigLength2):
	f = open('minimoFilter.sh', 'w')
	job=0
	for (key, vdict) in seeds.items():
		f.write(dirscr+'faLenFilter.py '+wd+'/fastq/'+key+'_contig2-contigs.fa '+wd+'/fastq/'+key+'_contig4 '+str(contigLength2)+' '+base+'_'+key+'\n')
	f.close()

def cap3MergeSinglet(contigLength2):
	f = open('cap3MergeSinglet.sh', 'w')
	job=0
	for (key, vdict) in seeds.items():
		f.write('cat '+wd+'/fastq/'+key+'_contig2.cap.singlets '+wd+'/fastq/'+key+'_contig2.cap.contigs > '+ wd+'/fastq/'+key+'_contig3\n')
		f.write(dirscr+'faLenFilter.py '+wd+'/fastq/'+key+'_contig3 '+ wd+'/fastq/'+key+'_contig4 '+str(contigLength2)+' '+base+'_'+key+'\n')
	f.close()

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
	if 'a' in assembly_para:
		r7=abyss_partition()
		sf.write('source abyss_partition.sh\nwait\n')
		sf.write('source abyss_combine.sh\nwait\n')
	moveContig1(assembly_para, r1,r2,r3,r4, r5,r6,r7,r8,r9,r10,r11)
	sf.write('source moveContigIndividual.sh\n')
	contig_out = wd+'/contig_individual/'
	sf.close()

def Assembly_combine(assembly_para):
	contig_out = wd+'/contig_'+assembly_para+'/'
	try: os.mkdir(contig_out)
	except: pass
	sf=open(contig_out+'assembly_combine.sh', 'w')
	combineContig2(contigLength1,assembly_para)
	os.system('cp combineContig2.sh '+contig_out+'\nwait\n')
	sf.write('source '+contig_out+'combineContig2.sh \nwait\n')
	if 'C' in assembly_para: 
		cap3()
		cap3MergeSinglet(contigLength2)
		os.system('cp cap3.sh '+contig_out+'\nwait\n')
		os.system('cp cap3MergeSinglet.sh '+contig_out+'\nwait\n')
		sf.write('source '+contig_out+'cap3.sh\nwait\n')
		sf.write('source '+contig_out+'cap3MergeSinglet.sh \nwait\n')
	elif 'O' in assembly_para:
		minimo()
		minimoFilter(contigLength2)
		os.system('cp minimo.sh '+contig_out+'\nwait\n')
		os.system('cp minimoFilter.sh '+contig_out+'\nwait\n')
		sf.write('source '+contig_out+'minimo.sh\nwait\n')
		sf.write('source '+contig_out+'minimoFilter.sh\nwait\n')

	moveContig2(assembly_para)
	sf.write('source '+contig_out+'moveContig.sh\n')
	sf.close()

def AddPipe(para, sf):
	assembly_para=para
	Assembly_combine(assembly_para)
	sf.write('source '+wd+'/contig_'+assembly_para+'/'+'assembly_combine.sh\n')

if __name__ == "__main__":
	configfile = sys.argv[1]
	readConfig(configfile)
	wd = os.path.abspath('.')
	base = os.path.basename(wd)
	try: 
		os.system('rm -rf '+wd+'/*.sh\n')
		os.system('rm -rf '+wd+' */\n')
	except: pass

	print 'path',wd
	print 'thread', thread
	print 'base', base
	print 'seeds', seeds
	print 'pairend', pair
	print 'metakmer', metakmer
	print 'soapkmer', soapkmer
	print 'abysskmer', abysskmer
	print 'LEN_DBG', contigLength1 #length before cap3
	print 'LEN_OLC', contigLength2 #length after cap3 before blast
	print 'mode', mode

	prep_reads(pair)
	clean_dir()

	if mode=='quick': Assembly_individual('Aa')
	else:Assembly_individual('SAVa')
	sf=open('ensemble.sh', 'w')
	sf.write('source prep_reads.sh\nwait\n')
	sf.write('source assembly_individual.sh\n')
	if mode == 'quick': AddPipe('AaO', sf)
	else: AddPipe('SAVaC', sf)
	sf.write('cat contig*/contig.log > contig.log\n')
	sf.close()
	try:  os.system('rm -rf '+wd+'/fastq/\n')
	except: pass