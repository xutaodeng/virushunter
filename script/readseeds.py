#!/usr/bin/env python
from collections import defaultdict
import os.path
import os
import string
import sys

def readSeeds(seedfile):
	seeds={}
	f=open(seedfile, 'r')
	data = f.readlines()
	data = ''.join(data).split('\"')
	print len(data)
	for i in xrange(len(data)/2):
		oldkey = '_'.join(data[2*i].strip().split('\t')[0].split())
		key=[]
		for letter in oldkey:
			if not (letter in string.letters or letter in string.digits):
				key.append('_')
			else:
				key.append(letter)
		key=''.join(key)
		value = data[2*i+1].split('>')
		j=0
		print i, key
		vdict ={}
		for j in xrange(1, len(value)):
			v = value[j]
			print '\t', j, 
			j+=1

			parts = v.split('\n', 1)
			id = '_'.join(parts[0].split())
			try: seq = ''.join(parts[1].split())
			except: seq=''
			print id, len(seq), seq
			vdict[id]=seq
		seeds[key]=vdict
	return seeds

def readTag(seedfile):
	seeds={}
	f=open(seedfile, 'r')
	for line in f:
		key, tag =line.strip().split()
		seeds[key]=tag
	return seeds
	

# now generate commandline sh for SOAPDENOVO
def soap():
	f = open('soap.sh', 'w')
	try: os.mkdir(os.path.abspath('.')+'/soap_config/')
	except: pass
	try: os.mkdir(os.path.abspath('.')+'/soap_out/')
	except: pass
	for (key, vdict) in seeds.items():
		ff = open(os.path.abspath(wd)+'/soap_config/'+key+'_soap.config', 'w')
		f.write('/mnt/cluster/tools/SOAPdenovo31mer all -K 31')
		f.write(' -s '+ os.path.abspath(wd)+'/soap_config/'+key+'_soap.config ')
		f.write(' -o '+ os.path.abspath(wd)+'/soap_out/'+key+'_soap\n')
		ff.write('max_rd_len=80\n[LIB]\n')
		ff.write('avg_ins=200\nreverse_seq=0\nasm_flags=3\nrank=1\n')
		ff.write('q1=fastq/'+key+'_1_sequence.txt\n'+'q2=fastq/'+key+'_2_sequence.txt\n')
		ff.close()
	f.close()

def splitQuery(n): #split fasta files into smaller chunks before virus blast
	f = open(os.path.abspath(wd)+'/splitQuery.sh', 'w')
	dbname= os.path.abspath('/mnt/cluster/xdeng/blastdb/virus/virus_prot') #virus blastdb
	for (key, vdict) in seeds.items():
		queryname = os.path.abspath(wd)+'/soap_out/'+key+'_soap.contig2'
		f.write('/mnt/cluster/xdeng/script/splitQuery.py '+queryname+' '+str(n)+'\n')
	f.close()

def mergeVirusXML(n):#merge blast XML files
	l=50
	f = open(os.path.abspath(wd)+'/mergeXML.sh', 'w')
	for (key, vdict) in seeds.items():
		outname = os.path.abspath(wd)+'/blast_virus_out/'+key+'_blast.xml'
		f.write('/mnt/cluster/xdeng/script/mergeXML.py '+outname+' '+str(n)+'\n')
		f.write('cat ')
		for i in xrange(n):
			queryname = os.path.abspath(wd)+'/soap_out/'+key+'_soap.contig2'
			parts=queryname.rsplit('.')
			subquery=parts[0]+'_'+str(i)+'.'+parts[1]
			subfilter=subquery.rsplit('.', 1)[0]+str(l)+'_filter.fa'
			f.write(subfilter+' ')
		f.write('> '+queryname.rsplit('.', 1)[0]+str(l)+'_filter.fa\n')
	f.close()

def blastVirusSplit(n):
	f = open(os.path.abspath(wd)+'/blast_virus.sh', 'w')
	dbname= os.path.abspath('/mnt/cluster/xdeng/blastdb/virus/virus_prot') #virus blastdb
	for (key, vdict) in seeds.items():
		queryname = os.path.abspath(wd)+'/soap_out/'+key+'_soap.contig2'
		outname = os.path.abspath(wd)+'/blast_virus_out/'+key+'_blast.xml'
		try: os.mkdir(os.path.abspath(wd)+'/blast_virus_out/')
		except: pass
		for i in xrange(n):
			f.write('/mnt/cluster/xdeng/script/blast_wrapper.py'+' blastx -num_threads 7 -evalue 0.01 -length 50 -outfmt 5 -dbsize 10000000 -searchsp 1000000000 -max_target_seqs 2')
			parts=queryname.rsplit('.')
			subquery=parts[0]+'_'+str(i)+'.'+parts[1]
			parts=outname.rsplit('.')
			subxml=parts[0]+'_'+str(i)+'.'+parts[1]
			f.write(' -db '+dbname+ ' -query '+ subquery+' -out '+subxml+'\n')
		#/mnt/cluster/xdeng/Eric/454/colorado/1.TCA.454Reads_W2_assembly/1.TCA.454Reads_W2_d_results/1.TCA.454Reads_W2_out.unpadded.fasta
	f.close()
	
def blastVirus():
	f = open(os.path.abspath(wd)+'/blast_virus.sh', 'w')
	dbname= os.path.abspath('/mnt/cluster/xdeng/blastdb/virus/virus_prot') #virus blastdb
	for (key, vdict) in seeds.items():
		queryname = os.path.abspath(wd)+'/soap_out/'+key+'_soap.contig2'
		outname = os.path.abspath(wd)+'/blast_virus_out/'+key+'_blast.xml'
		try: os.mkdir(os.path.abspath(wd)+'/blast_virus_out/')
		except: pass
		f.write('/mnt/cluster/xdeng/script/blast_wrapper.py'+' blastx -num_threads 7 -evalue 0.01 -length 50 -outfmt 5 -dbsize 10000000 -searchsp 1000000000 -max_target_seqs 2')
		f.write(' -db '+dbname+ ' -query '+ queryname+' -out '+outname+'\n')
		#/mnt/cluster/xdeng/Eric/454/colorado/1.TCA.454Reads_W2_assembly/1.TCA.454Reads_W2_d_results/1.TCA.454Reads_W2_out.unpadded.fasta
	f.close()

def blastVirusParser():
# now generate commandline sh for blast soap vs. virus, blastx
	f = open(os.path.abspath(wd)+'/blast_virus_parser.sh', 'w')
	dbname= os.path.abspath('/mnt/cluster/xdeng/blastdb/virus/virus_prot') #virus blastdb
	for (key, vdict) in seeds.items():
		queryname = os.path.abspath(wd)+'/soap_out/'+key+'_soap.contig2'
		faname=os.path.dirname(queryname)+'/'+os.path.basename(queryname).rsplit('.', 1)[0]+'50_filter.fa'
		xmlname = os.path.abspath(wd)+'/blast_virus_out/'+key+'_blast.xml'
		sigfaname =  os.path.dirname(queryname)+'/'+os.path.basename(queryname).rsplit('.', 1)[0]+'50_filterSig.fa' #output significant virus fa file for NR
		outname = os.path.abspath(wd)+'/blast_virus_out/'+key+'_blast.txt'
		f.write('/mnt/cluster/xdeng/script/blast_parser.py '+xmlname+ ' '+ faname+ ' '+ sigfaname+' > '+outname+'\n')
		#/mnt/cluster/xdeng/Eric/454/colorado/1.TCA.454Reads_W2_assembly/1.TCA.454Reads_W2_d_results/1.TCA.454Reads_W2_out.unpadded.fasta
	f.close()

def blastNR(nserver):
	f = open(os.path.abspath(wd)+'/blast_nr.sh', 'w')
	dbname= os.path.abspath('/mnt/cluster/xdeng/blastdb/nr/nr') #virus blastdb
	counter=0
	for (key, vdict) in seeds.items():
		serverTag = 'ssh bsidna'+str(counter%nserver+1)+' '
		queryname = os.path.abspath(wd)+'/soap_out/'+key+'_soap.contig2'
		sigfaname =  os.path.dirname(queryname)+'/'+os.path.basename(queryname).rsplit('.', 1)[0]+'50_filterSig.fa' #output significant virus fa file for NR
		outname = os.path.abspath(wd)+'/blast_nr_out/'+key+'_blast.xml'
		try: os.mkdir(os.path.abspath(wd)+'/blast_nr_out/')
		except: pass
		f.write(serverTag+' /mnt/cluster/xdeng/script/blast_wrapper.py '+' blastx -num_threads 8 -evalue 0.01 -outfmt 5 -dbsize 10000000 -searchsp 1000000000 -max_target_seqs 2')
		f.write(' -db '+dbname+ ' -query '+ sigfaname+' -out '+outname+' &\n')
		counter+=1
	f.close()

def blastNRFilter():
	f = open(os.path.abspath(wd)+'/blast_nr_filter.sh', 'w')
	f2 = open(os.path.abspath(wd)+'/blast_output_sort.sh', 'w')
	for (key, vdict) in seeds.items():
		queryname = os.path.abspath(wd)+'/soap_out/'+key+'_soap.contig2'
		faname=os.path.dirname(queryname)+'/'+os.path.basename(queryname).rsplit('.', 1)[0]+'50_filter.fa'
		nrXML = os.path.abspath(wd)+'/blast_nr_out/'+key+'_blast.xml'
		virusXML = os.path.abspath(wd)+'/blast_virus_out/'+key+'_blast.xml'
		try: os.mkdir(os.path.abspath(wd)+'/blast_filter_out/')
		except: pass
		outname = os.path.abspath(wd)+'/blast_filter_out/'+key+'_blast_filter.txt'
		outsorted = os.path.abspath(wd)+'/blast_filter_out/sorted_'+key+'_blast_filter.txt'
		outfa = os.path.abspath(wd)+'/blast_filter_out/sorted_'+key+'_blast_filter.fa'
		f.write('/mnt/cluster/xdeng/script/blast_filter_NR.py '+ virusXML + '  '+ nrXML+'  '+ faname+ '> '+outname+'\n')
		f2.write('/mnt/cluster/xdeng/script/blast_output_sort.py '+ outname +' '+outsorted+' '+outfa+'\n')
		#/mnt/cluster/xdeng/Eric/454/colorado/1.TCA.454Reads_W2_assembly/1.TCA.454Reads_W2_d_results/1.TCA.454Reads_W2_out.unpadded.fasta
	f.close()
	f2.close()

def combineContig_reads():
	f = open(os.path.abspath(wd)+'/combineContig_reads.sh', 'w')
	for (key, vdict) in seeds.items():
		f.write('cat '+os.path.abspath(wd)+'/fastq/'+key+'_1_sequence.txt ')
		f.write(os.path.abspath(wd)+'/fastq/'+key+'_2_sequence.txt ')
		f.write(' > '+os.path.abspath(wd)+'/fastq/'+key+'.fq \n')
		f.write('/mnt/cluster/xdeng/script/fq2fa.py ')
		f.write(os.path.abspath(wd)+'/fastq/'+key+'.fq '+os.path.abspath(wd)+'/fastq/'+key+'.fa 50 \n')
		queryname = os.path.abspath(wd)+'/soap_out/'+key+'_soap.contig'
		contigname = os.path.abspath(wd)+'/soap_out/'+key+'_soap.contig2' #contig combined with reads
		#faname=os.path.dirname(queryname)+'/'+os.path.basename(queryname).rsplit('.', 1)[0]+'50_filter.fa'
		#faout=os.path.dirname(queryname)+'/'+os.path.basename(queryname).rsplit('.', 1)[0]+'50_filter_out.fa'
		f.write('cat '+queryname+' '+os.path.abspath(wd)+'/fastq/'+key+'.fa > '+contigname +'\n')
		
	f.close()

if __name__ == "__main__":
	# tagFile=sys.argv[1] #'barcode.txt'
	# seeds=readTag(tagFile)
	try: 
		#seq1File = sys.argv[1] #'myseq.fastq'
		#seq2File = sys.argv[2] #'myseq.fastq'
		tagFile='barcode.txt'
		seeds=readTag(tagFile)
		#seeds=readSeeds(tagFile)
	except:
		seeds = {'myseq':'myseq'}
	
	wd = os.path.dirname('.')
	print wd

	n=5
	soap()
	combineContig_reads()
	splitQuery(n)
	#blastVirus()
	blastVirusSplit(n)
	mergeVirusXML(n)
	blastVirusParser()
	blastNR(30)
	blastNRFilter()

	# os.system('PriceTI  -fp s1.fq s2.fq 300 -icf seeds.fa 1 1 5  -nc 30 -dbmax 72 -mol 35 -tol 20 -mpi 80 -target 90 2 2 2 -a 7 -o price.fa')
	# os.system('cat read1_sequence.txt read2_sequence.txt seq.fq')
	# os.system('clone_reads_rm.py seq.fq clone.fastq')
	# os.system('trim.py clone.fastq trim.fastq')
	# os.system('norm_fq_pair.py trim.fastq seq1.fq seq2.fq')
	# os.system('fq_pair_addIndex.py index1.txt index2.txt seq1.fq seq2.fq s1.fq s2.fq')
	# os.system('debarcode_pair.py s1.fq s2.fq barcode.txt')
	# os.system('sh '+os.path.abspath(wd)+'/soap.sh')
	# os.system('sh '+os.path.abspath(wd)+'/combineContig_reads.sh')
	# os.system('sh '+os.path.abspath(wd)+'/blast_virus.sh')
	# os.system('sh '+os.path.abspath(wd)+'/blast_virus_parser.sh')
	# os.system('sh '+os.path.abspath(wd)+'/blast_nr.sh')
	# os.system('sh '+os.path.abspath(wd)+'/blast_nr_filter.sh')
	# os.system('sh '+os.path.abspath(wd)+'/blast_output_sort.sh')
	
