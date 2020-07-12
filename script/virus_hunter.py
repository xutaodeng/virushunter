#!/usr/bin/env python
from collections import defaultdict
import os.path
import os
import string
import sys
import itertools
import subprocess

# samtools lcurses to lncurses
# find . -iname "*.gz" -type f -exec /bin/mv {} /mnt/cluster2/xdeng/SatishRNASeq/ \;
# https://stat.ethz.ch/pipermail/r-help/2006-April/103372.html

# rsync -ru /source/directory/* username@domain.net:/destination/directory
##################################################################################################
# 1. prepare human and bacteria filters for bowtie
##################################################################################################
# wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz

# wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.zip
# wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_prot.zip
# wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxcat.zip
# wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip
# 7za e taxcat.zip
# 7za e gi_taxid_nucl.zip
# 7za e gi_taxid_prot.zip
# 7za e taxdmp.zip

## bacteria preparation
# /mnt/cluster/xdeng/nt/nt_extract_bac.py

# ssh bsidna3 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria1.fa /mnt/cluster/xdeng/nt/Bacteria1 &
# ssh bsidna4 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria2.fa /mnt/cluster/xdeng/nt/Bacteria2 &
# ssh bsidna5 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria3.fa /mnt/cluster/xdeng/nt/Bacteria3 &
# ssh bsidna6 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria4.fa /mnt/cluster/xdeng/nt/Bacteria4 &
# ssh bsidna7 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria5.fa /mnt/cluster/xdeng/nt/Bacteria5 &
# ssh bsidna8 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria6.fa /mnt/cluster/xdeng/nt/Bacteria6 &
# ssh bsidna10 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria7.fa /mnt/cluster/xdeng/nt/Bacteria7 &
# ssh bsidna11 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria8.fa /mnt/cluster/xdeng/nt/Bacteria8 &
# ssh bsidna14 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria9.fa /mnt/cluster/xdeng/nt/Bacteria9 &
# ssh bsidna15 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria10.fa /mnt/cluster/xdeng/nt/Bacteria10 &
#wait
# ssh bsidna3 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria1.fa /mnt/cluster/xdeng/nt/Bacteria11 &
# ssh bsidna4 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria2.fa /mnt/cluster/xdeng/nt/Bacteria12 &
# ssh bsidna5 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria3.fa /mnt/cluster/xdeng/nt/Bacteria13 &
# ssh bsidna6 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria4.fa /mnt/cluster/xdeng/nt/Bacteria14 &
# ssh bsidna7 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria5.fa /mnt/cluster/xdeng/nt/Bacteria15 &
# ssh bsidna8 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria6.fa /mnt/cluster/xdeng/nt/Bacteria16 &
# ssh bsidna10 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria7.fa /mnt/cluster/xdeng/nt/Bacteria17 &
# ssh bsidna11 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria8.fa /mnt/cluster/xdeng/nt/Bacteria18 &
# ssh bsidna14 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria9.fa /mnt/cluster/xdeng/nt/Bacteria19 &
# ssh bsidna15 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria10.fa /mnt/cluster/xdeng/nt/Bacteria20 &
# ssh bsidna19 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria1.fa /mnt/cluster/xdeng/nt/Bacteria21 &
# ssh bsidna20 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria2.fa /mnt/cluster/xdeng/nt/Bacteria22 &
# ssh bsidna23 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria3.fa /mnt/cluster/xdeng/nt/Bacteria23 &
# ssh bsidna24 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria4.fa /mnt/cluster/xdeng/nt/Bacteria24 &
# ssh bsidna25 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria5.fa /mnt/cluster/xdeng/nt/Bacteria25 &
# ssh bsidna26 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria6.fa /mnt/cluster/xdeng/nt/Bacteria26 &
# ssh bsidna30 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria7.fa /mnt/cluster/xdeng/nt/Bacteria27 &



# servers=[  'bsidna4', \
 # 'bsidna5', 'bsidna6','bsidna7','bsidna8', \
  # 'bsidna9', 'bsidna10','bsidna11', \
  # 'bsidna14', 'bsidna15', 'bsidna19', 'bsidna20', \
  # 'bsidna23', 'bsidna24', 'bsidna25', 'bsidna26', \
# 'bsidna27','bsidna28','bsidna29', 'bsidna30'


# /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/insect/insect.fasta /mnt/cluster/xdeng/insect/insect
# /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/zebrafish/zebra.fasta /mnt/cluster/xdeng/zebrafish/zebra
# /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/mosquito/mosquito.fasta /mnt/cluster/xdeng/mosquito/mosquito
#wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
#wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/mrna.fa.gz
# cd hg38
# zcat hg38.fa.gz mrna.fa.gz > mrnadna.fa
# ssh bsidna12 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2 /mnt/cluster/xdeng/hg38/mrnadna.fa /mnt/cluster/xdeng/hg38/mrnadna_bowtie &


# #0.	Add 7zip to path:
# #a.	/mnt/cluster/tools/p7zip_9.20.1/bin

# #1.	Back up the old nr and virus database:

# #mv /mnt/san/cluster/xdeng/blastdb/nr /mnt/san/cluster/xdeng/blastdb/nr_today
# #mv /mnt/san/cluster/xdeng/blastdb/virus /mnt/san/cluster/xdeng/blastdb/virus_today

# #2.	Download the virus genome and nr to current directory /mnt/san/cluster/xdeng/blastdb/

# wget ftp://ftp.ncbi.nih.gov/refseq/release/viral/viral.1.protein.faa.gz
# wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
## wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/*/latest_assembly_versions/*/*_genomic.fna.gz

# cd refseq_download
# wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/*protein.faa.gz -P ./refseq_download
# zcat ./refseq_download/*.gz > refseq.fa

# #3.	Update taxonomy:
# wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_prot.zip
# wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxcat.zip
# 7za e taxcat.zip
# 7za e gi_taxid_prot.zip
# 7za e nr.zip
# #4.	Create new empty directory
# mkdir /mnt/san/cluster/xdeng/blastdb/nr
# mkdir /mnt/san/cluster/xdeng/blastdb/virus
# #5.	Process the downloaded files 
                # python nr_virus.py
# #6.	Run makeblastdb
# segmasker -locut 1.9 -hicut 2.5 -in virus.fa -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out virus_mask.asnb
# makeblastdb -in virus.fa -dbtype prot -parse_seqids -mask_data virus_mask.asnb -out virus_mask
# /mnt/cluster/xdeng/tools/diamond/diamond makedb --in virus.fa -d virusdiamond
# /mnt/cluster/xdeng/tools/diamond/diamond makedb --in diamond.fa -d diamond

#/mnt/cluster/xdeng/taxon> cp viral_pasteur.fa ../blastdb/
# segmasker -locut 1.9 -hicut 2.5 -in viral_pasteur.fa -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out viral_pasteur_mask.asnb
# makeblastdb -in viral_pasteur.fa -dbtype prot -parse_seqids -mask_data viral_pasteur_mask.asnb -out viral_pasteur_mask


 # segmasker -locut 1.9 -hicut 2.5 -in phan_phage.fa -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out phan_phage_mask.asnb
 # makeblastdb -in phan_phage.fa -dbtype prot -parse_seqids -mask_data phan_phage_mask.asnb -out phan_phage_mask



# dustmasker -in Pool31to35_m.fasta -infmt fasta -outfmt fasta -out test_out.fasta

# segmasker -locut 1.9 -hicut 2.5 -in phageHERV.fa -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out phageHERV_mask.asnb
# makeblastdb -in phageHERV.fa -dbtype prot -parse_seqids -mask_data phageHERV_mask.asnb -out phageHERV_mask
# #segmasker -locut 1.9 -hicut 2.5 -in nvrefseq.fa -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out nvrefseq_mask.asnb
# #makeblastdb -in nvrefseq.fa -dbtype prot -parse_seqids -mask_data nvrefseq_mask.asnb  -out nvrefseq_mask
# makeblastdb -in nvrefseq.fa -dbtype prot -parse_seqids -out nvrefseq

 # prerapsearch -d virus.fa -n virusrap
 # prerapsearch -d nvnr.fa -n nvrap
 #rapsearch -q 4440037.3.dna.fa -d nogCOGdomN95_db -o 4440037.3.dna-vs-nogCOGdomN95 -z 8 -b 1 -v 1 -x t

# #7.	Change mode
# #a.	sudo chmod 777 * -R

 # md5sum -b *.gz > md5.txt
 # tar cvzf sra_metagenomics.tgz *
 # ftp ftp-private.ncbi.nih.gov
# Address: ftp-private.ncbi.nih.gov
# Login: sra
# Password: 03o4!srRy!
 # put sra_metagenomics.tgz
 # bye
# ======================================================================================================
#wget -r -l1 --no-parent -A.gz  --user="eric@oucru.org" --password="Oucru@@2018share" ftp://files.oucru.org//ERIC/Miseq/34th/ 
#Instruction for running virus discovery pipeline
# 0 (set up one time only):
                # setup path for the following tools
                # /mnt/cluster/xdeng/script
                # /usr/bin/
                # /mnt/cluster2/xdeng/tools
# 1. download the data into directory:
                # /mnt/cluster/xdeng/Eric/NewProject/fastq
        # wget -r -l1 --no-parent -A.gz --user='chiulab' --password='!wantdata' vddc.ucsf.edu/~miseq/miseq_data/131021_LL/
		# wget -r -l1 --no-parent --user='xutao' --password='den0v0' 	chiulab.ucsf.edu/~chiulab/denovo-assembly/additional_datasets_denovo_assembly_5-28-14.tgz


                # change the mode to be rwe for all files in the project directory
# 2. create sample file inside fastq folder:
                # ls -1 *.gz >samples.txt
# 3. go to project directory: 
                # cd ..
                # run readseeds2.py to generate pipeline script
                # sudo chmod 777 * -R
                  # #sh pipeline_run.sh

# Also reinstall blast software every Jan and July or if there is a new version
# ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/



RAM='25G'


servers=[  'bsidna4', \
 'bsidna5', 'bsidna6','bsidna7','bsidna8', \
  'bsidna9', 'bsidna10','bsidna11', \
  'bsidna14', 'bsidna15', 'bsidna19', 'bsidna20', \
  'bsidna23', 'bsidna24', 'bsidna25', 'bsidna26', \
'bsidna27','bsidna28','bsidna29', 'bsidna30'
]
#servers=[  'bsidna35']

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

# servers=['bsidna3']
trinitypath='/mnt/cluster/xdeng/tools/trinityrnaseq-2.2.0/'
#spadepath='/mnt/cluster/xdeng/tools/SPAdes-3.8.0-Linux/bin/'
spadepath='/mnt/cluster/tools/SPAdes-3.11.1-Linux/bin/'
#spadepath='/mnt/cluster/tools/SPAdes-3.13.1-Linux/bin/'
#diamondpath='/mnt/cluster/tools/diamond'
cap3path = '/mnt/cluster/xdeng/tools/CAP3/cap3'
soappath='/mnt/cluster/xdeng/tools/SOAPdenovo2-src-r240/SOAPdenovo-63mer'
velvetg='/mnt/cluster/tools/MetaVelvet-1.2.02/velvetg'
velveth='/mnt/cluster/tools/MetaVelvet-1.2.02/velveth'
abysspath='/mnt/cluster/tools/abyss/bin/abyss-pe'
meta_velvetg="/mnt/cluster/tools/MetaVelvet-1.2.02/meta-velvetg"
bowtiepath='/mnt/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2' 
blastxpath='/mnt/cluster/xdeng/tools/ncbi-blast-2.2.31+/bin/blastx'
blastnpath='/mnt/cluster/xdeng/tools/ncbi-blast-2.2.31+/bin/blastn'
dustmaskerpath='/mnt/cluster/xdeng/tools/ncbi-blast-2.2.31+/bin/dustmasker'
samtools='/mnt/cluster/xdeng/tools/samtools-1.2/samtools'
rapsearch='/mnt/cluster/xdeng/tools/RAPSearch2.23_64bits/bin/rapsearch'
prerapsearch='/mnt/cluster/xdeng/tools/RAPSearch2.23_64bits/bin/prerapsearch'
newdiamondpath='/mnt/cluster/xdeng/tools/diamond/diamond' #new version of diamond can generate xml blast format
diamondpath=newdiamondpath

diamonddb='/mnt/cluster/xdeng/blastdb/diamond.dmnd'
bowtieindexpath='/mnt/cluster/xdeng/hg38/mrnadna_bowtie'
			
bowtiebacs= [
			'/mnt/cluster/xdeng/nt/Bacteria1', \
			'/mnt/cluster/xdeng/nt/Bacteria2', \
			'/mnt/cluster/xdeng/nt/Bacteria3', \
			'/mnt/cluster/xdeng/nt/Bacteria4', \
			'/mnt/cluster/xdeng/nt/Bacteria5', \
			'/mnt/cluster/xdeng/nt/Bacteria6', \
			'/mnt/cluster/xdeng/nt/Bacteria7', \
			'/mnt/cluster/xdeng/nt/Bacteria8', \
			'/mnt/cluster/xdeng/nt/Bacteria9', \
			'/mnt/cluster/xdeng/nt/Bacteria10', \
			'/mnt/cluster/xdeng/nt/Bacteria11', \
			'/mnt/cluster/xdeng/nt/Bacteria12', \
			'/mnt/cluster/xdeng/nt/Bacteria13', \
			'/mnt/cluster/xdeng/nt/Bacteria14', \
			'/mnt/cluster/xdeng/nt/Bacteria15', \
			'/mnt/cluster/xdeng/nt/Bacteria16', \
			'/mnt/cluster/xdeng/nt/Bacteria17', \
			'/mnt/cluster/xdeng/nt/Bacteria18', \
			'/mnt/cluster/xdeng/nt/Bacteria19', \
			'/mnt/cluster/xdeng/nt/Bacteria20', \
			'/mnt/cluster/xdeng/nt/Bacteria21', \
			'/mnt/cluster/xdeng/nt/Bacteria22', \
			'/mnt/cluster/xdeng/nt/Bacteria23', \
			'/mnt/cluster/xdeng/nt/Bacteria24', \
			'/mnt/cluster/xdeng/nt/Bacteria25', \
			'/mnt/cluster/xdeng/nt/Bacteria26', \
			'/mnt/cluster/xdeng/nt/Bacteria27'
			
			#'/mnt/cluster/xdeng/malaria/malaria'
			#'/mnt/cluster/xdeng/insect/insect',
			#'/mnt/cluster/xdeng/zebrafish/zebra' #,\
			#'/mnt/cluster/xdeng/mosquito/mosquito'
			]
bowtieNTIndex= ['/mnt/cluster/xdeng/taxon/nt_1', \
				'/mnt/cluster/xdeng/taxon/nt_2', \
				'/mnt/cluster/xdeng/taxon/nt_3', \
				'/mnt/cluster/xdeng/taxon/nt_4', \
				'/mnt/cluster/xdeng/taxon/nt_5', \
				'/mnt/cluster/xdeng/taxon/nt_6', \
				'/mnt/cluster/xdeng/taxon/nt_7', \
				'/mnt/cluster/xdeng/taxon/nt_8', \
				'/mnt/cluster/xdeng/taxon/nt_9', \
				'/mnt/cluster/xdeng/taxon/nt_10', \
				'/mnt/cluster/xdeng/taxon/nt_11', \
				'/mnt/cluster/xdeng/taxon/nt_12', \
				'/mnt/cluster/xdeng/taxon/nt_13', \
				'/mnt/cluster/xdeng/taxon/nt_14'
			]
#including 3 insect genomes
diamondvirusdb='/mnt/cluster/xdeng/blastdb/virusdiamond.dmnd'
hmmannot='/mnt/cluster/xdeng/blastdb/vFam/annot2014.txt'
virusdbpath='/mnt/cluster/xdeng/blastdb/virus_mask'
virusdbpathDNA='/mnt/cluster/xdeng/blastdb/virus_DNA_mask' 
#virusdbpath='/mnt/cluster/xdeng/blastdb/viral_pasteur_mask'
phagedbpath='/mnt/cluster/xdeng/blastdb/phage_mask'
virusphagedbpath='/mnt/cluster/xdeng/blastdb/virusphage_mask'
#phan_phagedbpath='/mnt/cluster/xdeng/blastdb/phan_phage_mask'
#nvnrdbpath='/mnt/cluster/xdeng/blastdb/nvrefseq'
nvnrdbpath='/mnt/cluster/xdeng/blastdb/nvnr'
vfampath='/mnt/cluster/xdeng/blastdb/vFam/vFam-A_2014.hmm'
rapdbnr='/mnt/cluster/xdeng/blastdb/nvrap'
nvnrfa = '/mnt/cluster/xdeng/blastdb/nvrefseq.fa'
virusfa = '/mnt/cluster/xdeng/blastdb/virus.fa'
dirscr='/mnt/cluster/xdeng/script/'
Raytool='/mnt/cluster2/xdeng/tools/Ray2.3/Ray'
picard = '/mnt/cluster/xdeng/tools/picard-tools-1.136/'
nservers=len(servers)
nsuperservers=len(superservers)

#hmmerpath = '/mnt/cluster/xdeng/tools/hmmer-3.1b2-linux-intel-x86_64/binaries/phmmer'
hmmerpath = '/mnt/cluster/xdeng/tools/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmsearch'
# -o  hello  tutorial/HBB_HUMAN tutorial/globins45.fa


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

def Clark():
	global seeds, pair
	try: os.mkdir('clark_out')
	except: pass
	try: os.mkdir(wd+'/'+base+'/clark')
	except: pass
	try: os.mkdir(wd+'/'+base+'/clark/fasta')
	except: pass
	f = open('clark.sh', 'w')
	f2 = open('clark_process.sh', 'w')
	job=0
	for (key, fqfiles) in seeds.items():
		serverTag = 'ssh '+superservers[job%nsuperservers]
		serverTag2 = 'ssh '+servers[job%nservers]
		# fq = [wd+'/fastq/'+fqfile for fqfile in fqfiles]
		# for fqq in fq:
		f.write(serverTag+' /mnt/cluster/tools/CLARKSCV1.2.3.2/exe/CLARK -k 20 -n 48 -T /mnt/cluster/xdeng/taxon/target1.txt -D /mnt/cluster/xdeng/taxon/CLARK_DB1/ -O '+wd+'/fastq/'+key+'.fa -R '+wd+'/clark_out/'+key+'_clark1 &\n')
		f2.write('cat '+wd+'/clark_out/'+key+'*.csv > ' + wd+'/clark_out/'+key+'.csvv \n')
		indexfile='/mnt/cluster/xdeng/taxon/species_index.txt'
		f2.write(dirscr+'clark_result.py '+indexfile+' '+wd+'/clark_out/'+key+'.csvv '+wd+'/'+base+'/clark/'+key+'.count '+wd+'/fastq/'+key+'.fa '+wd+'/'+base+'\n')
		f2.write(dirscr+'clark_html.py '+wd+'/'+base+'/clark/'+key+'.count '+wd+'/'+base+'/clark/'+key+'.count.csv '+wd+'/'+base+'/clark/'+key+'.html \n')
		job+=1
		if job%nsuperservers==0:  f.write('wait\n')
	key='all'
	#f2.write('cat '+wd+'/clark_out/*.csvv > ' + wd+'/clark_out/all.csvvv \n')
	#f2.write(dirscr+'clark_result.py '+indexfile+' '+wd+'/clark_out/all.csvvv '+wd+'/'+base+'/clark/all.count '+wd+'/fastq/all.fa '+wd+'/'+base+'\n')
	#f2.write(dirscr+'clark_html.py '+wd+'/'+base+'/clark/'+key+'.count '+wd+'/'+base+'/clark/'+key+'.count.csv '+wd+'/'+base+'/clark/'+key+'.html \n')

	for (key, fqfiles) in seeds.items():
		serverTag = 'ssh '+superservers[job%nsuperservers]
		# fq = [wd+'/fastq/'+fqfile for fqfile in fqfiles]
		# for fqq in fq:
		f.write(serverTag+' /mnt/cluster/tools/CLARKSCV1.2.3.2/exe/CLARK -k 20 -n 48 -T /mnt/cluster/xdeng/taxon/target2.txt -D /mnt/cluster/xdeng/taxon/CLARK_DB2/ -O '+wd+'/fastq/'+key+'.fa -R '+wd+'/clark_out/'+key+'_clark2 &\n')
		job+=1
		if job%nsuperservers==0:  f.write('wait\n') 
	for (key, fqfiles) in seeds.items():
		serverTag = 'ssh '+superservers[job%nsuperservers]
		# fq = [wd+'/fastq/'+fqfile for fqfile in fqfiles]
		# for fqq in fq:
		f.write(serverTag+' /mnt/cluster/tools/CLARKSCV1.2.3.2/exe/CLARK -k 20 -n 48 -T /mnt/cluster/xdeng/taxon/target3.txt -D /mnt/cluster/xdeng/taxon/CLARK_DB3/ -O '+wd+'/fastq/'+key+'.fa -R '+wd+'/clark_out/'+key+'_clark3 &\n')
		job+=1
		if job%nsuperservers==0:  f.write('wait\n') 

	# job=0
	# seeds2={}
	# for (key, fqfiles) in seeds.items():
		# serverTag = 'ssh '+servers[job%nservers]
		# f.write(serverTag+' cat '+wd+'/fastq/'+key+'.notCombined_1.fastq '+ wd+'/fastq/'+key+'.notCombined_2.fastq '+wd+'/fastq/'+key+'.extendedFrags.fastq > '+wd+'/fastq/'+key+'_flash.fq &\n')
		# seeds2[key]=[key+'_flash.fq']
		# job+=1
		# if job%nservers==0:  f.write('wait\n')
	# seeds=seeds2
	# for key in seeds:
		# print key, seeds[key]
	# pair=False
	# f.close()
	#return (seeds, pair)
	
# /mnt/cluster/tools/CLARKSCV1.2.3.2/exe/CLARK -k 20 -n 20 -T /mnt/cluster/xdeng/taxon/targets_addresses1.txt -D /mnt/cluster/xdeng/taxon/CLARK_DB1/ -O /mnt/cluster2/Vhunt/20170822_Anh_Miseq30th/fastq/4_S4_L001_R1_001.fastq.gz -R ./results1 

# /mnt/cluster/tools/CLARKSCV1.2.3.2/exe/CLARK -k 20 -n 20 -T /mnt/cluster/xdeng/taxon/targets_addresses2.txt -D /mnt/cluster/xdeng/taxon/CLARK_DB2/ -O /mnt/cluster2/Vhunt/20170822_Anh_Miseq30th/fastq/4_S4_L001_R1_001.fastq.gz -R ./results2 

# /mnt/cluster/tools/CLARKSCV1.2.3.2/exe/CLARK -k 20 -n 20 -T /mnt/cluster/xdeng/taxon/targets_addresses3.txt -D /mnt/cluster/xdeng/taxon/CLARK_DB3/ -O /mnt/cluster2/Vhunt/20170822_Anh_Miseq30th/fastq/4_S4_L001_R1_001.fastq.gz -R ./results 



def bam2sam(seedfile='bam/bams.txt'): #TCGA
	sams=[]
	try: os.mkdir('fastq')
	except: pass
	f=open(seedfile, 'r')
	of=open('bam.map', 'w')
	of2=open('bam2sam.sh', 'w')
	of3=open('samfilterfq.sh', 'w')
	i=0
	job=0
	for line in f:
		serverTag = 'ssh '+servers[job%nservers]
		i+=1
		sam=wd+'/bam/sam'+str(i)+'tmp.sam'
		fq=wd+'/fastq/sam'+str(i)+'.fastq'
		bam=wd+'/bam/'+line.strip()
		of2.write(serverTag+' '+samtools + ' view  -h -o '+sam+' '+bam+' &\n')
		of3.write(serverTag+' '+dirscr+'samfilterfq.py '+sam+' '+fq+' unmap&\n')
		of.write(bam+'\t'+sam+'\n')
		sams.append(sam)
		job+=1
		if job%nservers==0: 
			of2.write('wait\n')
			of3.write('wait\n')
	f.close()
	of.close()
	of2.close()
	of3.close()
	return sams

def tar2fastq(seedfile='fqtar/tars.txt'): #TCGA
	f=open(seedfile, 'r')
	of2=open('tar2fq.sh', 'w')
	job=0
	outdir=wd+'/fqtar/'
	for line in f:
		serverTag = 'ssh '+servers[job%nservers]
		tar=wd+'/fqtar/'+line.strip()
		of2.write(serverTag+' '+'tar zxvf '+tar+' -C '+outdir+' &\n')
		job+=1
		if job%nservers==0: 
			of2.write('wait\n')
	f.close()
	of2.close()

# def readSeeds2(seedfile='fastq/samples.txt'):
	# seeds=defaultdict(list)
	# f=open(seedfile, 'r')
	# for line in f:
		# key2=line.strip().rsplit('.', 1)[0]
		# key1=line.strip().rsplit('_',4)[0]
		# if len(key1)<len(key2):
			# key=key1
		# else:
			# key=key2
		# key=line.strip().rsplit('_',2)[0]
		# ss=line.strip()
		# if line.strip().endswith('fasta'): ss = ss.replace('fasta', 'fastq')
		# elif line.strip().endswith('sam'): ss = ss.replace('sam', 'fastq')
		# elif line.strip().endswith('bam'): ss = ss.replace('bam', 'fastq')
		# if len(seeds[key])<2: seeds[key].append(ss)
		# else:
			# key=key+'.2'
			# seeds[key].append(ss)
	# f.close()
	# return seeds


def mergeAB(reAssemble='no'):
	global seeds
	stats=defaultdict()
	seedfile='fastq/samples.txt'
	f=open(seedfile, 'r')
	seeds=defaultdict(list)
	of=open('merge.sh', 'w')
	for line in f:
		key2=line.strip().split('-', 1)[0]
		print key2
		key1=line.strip().rsplit('_',4)[0]
		if len(key1)<len(key2):
			key=key1
		else:
			key=key2
		#key=line.strip().rsplit('_',2)[0]
		ss=line.strip()
		k0=key
		if '_R1_' in line: key=key+'_R1'
		if '_R2_' in line: key=key+'_R2'

		seeds[key].append(ss)
		stats[key]=defaultdict(list)
	for key in seeds.keys():
		file=seeds[key][0]
		k0=file.strip().split('-', 1)[0]
		f2=file.strip().split('/', 1)[1]
		cmd='cat '+' '.join(seeds[key])+' > '+ k0+'_'+f2
		of.write( cmd+'\n')
	of.close()

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
	if reAssemble=='reAssemble': stats['reAssemble']=defaultdict(list)
	f.close()
	return stats, seeds
	
def readSeeds2(reAssemble='no'):
	global seeds
	stats=defaultdict()
	seedfile='fastq/samples.txt'
	f=open(seedfile, 'r')
	seeds=defaultdict(list)
	for line in f:
		key2=line.strip().rsplit('.', 1)[0]
		#key1=line.strip().rsplit('_',4)[0]
		key1=line.strip().rsplit('_',2)[0]
		print key1
		if len(key1)<len(key2):
			key=key1
		else:
			key=key2
		key=line.strip().split('.')[0].split('_')[1]
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
	if reAssemble=='reAssemble': stats['reAssemble']=defaultdict(list)
	f.close()
	return stats, seeds

def mergePairFq():
	global seeds, pair
	f = open('mergePairFq.sh', 'w')
	job=0
	for (key, fqfiles) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		fq = [wd+'/fastq/'+fqfile for fqfile in fqfiles]
		f.write(serverTag+' /mnt/cluster/xdeng/tools/FLASH-1.2.11/flash -M 250 -o '+key+' -d '+wd+'/fastq/ '+' '.join(fq)+' &\n')
		job+=1
		if job%nservers==0:  f.write('wait\n') 
	f.write('wait\n')
	job=0
	seeds2={}
	for (key, fqfiles) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		f.write(serverTag+' cat '+wd+'/fastq/'+key+'.notCombined_1.fastq '+ wd+'/fastq/'+key+'.notCombined_2.fastq '+wd+'/fastq/'+key+'.extendedFrags.fastq > '+wd+'/fastq/'+key+'_flash.fq &\n')
		seeds2[key]=[key+'_flash.fq']
		job+=1
		if job%nservers==0:  f.write('wait\n')
	seeds=seeds2
	for key in seeds:
		print key, seeds[key]
	pair=False
	f.close()
	#return (seeds, pair)

def skipbowtie():
	f = open('skipbowtie.sh', 'w')
	job=0
	for (key, fqfiles) in seeds.items():
		i=1
		for fqfile in fqfiles:
			serverTag = 'ssh '+servers[job%nservers]
			index=key+'_'+str(i)
			fqfil=index+'.fil'
			if fqfile.endswith('.gz'): f.write(serverTag+' zcat '+wd+'/fastq/'+fqfile +' > '+wd+'/fastq/'+fqfil+' &\n')
			else: f.write(serverTag+' cat '+wd+'/fastq/'+fqfile +' > '+wd+'/fastq/'+fqfil+' &\n')
			i+=1
			job+=1
			if job%nservers==0:  f.write('wait\n')
	if job%nservers==0:  f.write('wait\n')
	f.close()

def fa2fq():
	f = open('fa2fq.sh', 'w')
	for (key, fqfiles) in seeds.items():
		for fqfile in fqfiles:
			fafile = fqfile.replace('fastq', 'fasta')
			f.write(dirscr+'fa2fq2.py '+wd+'/fastq/'+fafile+' '+wd+'/fastq/'+fqfile+'\n') #quality trimming
	f.close()
	
def bam2fq():
	f = open('bam2fq.sh', 'w')
	job=0
	for (key, fqfiles) in seeds.items():
		for fqfile in fqfiles:
			serverTag = 'ssh '+servers[job%nservers]
			bamfile = fqfile.replace('fastq', 'bam')
			sam2fq = 'java  -Xms1024m -Xmx8g -jar  '+picard+'SamToFastq.jar INPUT='+wd+'/fastq/'+bamfile+' FASTQ='+wd+'/fastq/'+fqfile+' VALIDATION_STRINGENCY=SILENT'
			f.write(serverTag+' '+sam2fq+' &\n')
			if job%nservers==0:  f.write('wait\n')
			job+=1
	f.close()

def zip2fq():
	f = open('zip2fq.sh', 'w')
	for (key, fqfiles) in seeds.items():
		for fqfile in fqfiles:
			zip2fq = 'unzip '+wd+'/fastq/'+fqfile
			f.write(zip2fq+'\n')
	f.close()

def skip_adaptor():
	f = open('skipadaptor.sh', 'w')
	job=0
	for (key, fqfiles) in seeds.items():
		i=1
		for fqfile in fqfiles:
			serverTag = 'ssh '+servers[job%nservers]
			index=key+'_'+str(i)
			fid = wd+'/fastq/'+index
			if dedup: fqin = fid+'.dup'
			else: fqin=fid+'.fil'
			f.write(serverTag+' '+dirscr+'recodeID.py '+fqin+' '+fid+'.trim '+base+'_'+key+ ' '+str(i)+' &\n') #recode ID to line number
			i+=1
			job+=1
			if job%nservers==0:  f.write('wait\n')
	f.close()

def sampleFastq():
	f = open('sampleFastq.sh', 'w')
	for (key, fqfiles) in seeds.items():
		i=1
		fqs=[]
		for fqfile in fqfiles:
			fqs.append(fqfile)
			i+=1
		f.write(dirscr+'sampleFastq.py '+wd+'/fastq/'+fqs[0]+' '+wd+'/fastq/'+fqs[1]+' '+wd+'/fastq/out'+fqs[0][:-3]+' '+wd+'/fastq/out'+fqs[1][:-3]+' 75000\n') #quality trimming
	f.close()
	
def prepBlastFile():
	try: os.mkdir(wd+'/'+base+'/blast/')
	except: pass
	f = open('prepBlastFile.sh', 'w')
	#path='fastq/'
	#f.write('cd fastq/ \n')
	bpath=wd+'/'+base+'/blast/'
	directory = os.path.basename(wd) #all blast db
	job=0
	alldb=[]
	
	for (key, fqfiles) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		job+=1
		f.write(serverTag+' "cd '+wd+'/fastq/ &&')
		f.write('makeblastdb -dbtype nucl -parse_seqids -in '+key+'.fa -out '+key+'_blastdb "&\n')
		alldb.append(key+'_blastdb')
		if job%nservers==0:  f.write('wait\n')
	#f.write('blastdb_aliastool -dblist '+'\"'+' '.join(alldb)+'\" -dbtype nucl -out '+directory+'_blastdb -title \"'+directory+'_blastdb\"\n')
	f.write('wait\n')
	f.write('mv fastq/*_blastdb* '+bpath+'\n')
	#f.write('cd ..\n')
	f.write(dirscr+'blastAlias.py '+directory+' '+bpath+'\n')
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
			f.write('makeblastdb -dbtype nucl -parse_seqids -in '+index+'.fa -out '+index+'_blastdb && ')
			pairdb.append(index+'_blastdb')
			i+=1
		f.write('blastdb_aliastool -dblist '+'\''+' '.join(pairdb)+'\' -dbtype nucl -out '+key+'_blastdb -title \''+key+'_blastdb\' "&\n')
		if job%nservers==0:  f.write('wait\n')
		alldb.append(key+'_blastdb')
	f.write('wait\n')
	# f.write('blastdb_aliastool -dblist '+'\''+' '.join(alldb)+'\' -dbtype nucl -out '+directory+'_blastdb -title \''+directory+'_blastdb\' \n')
	f.write('mv fastq/*_blastdb* '+bpath+'\n')
	#f.write('cd ..\n')
	f.close()

def trim():
	f = open('clonetrim.sh', 'w')
	f1 = open('clonefq2fa.sh', 'w')
	f11=open('original_fq2fa.sh', 'w')
	f2 = open('blast_adaptor.sh', 'w')
	f3 = open('blasttrim.sh', 'w')
	f4 = open('qualitytrim.sh', 'w')
	f5 = open('fq_clean.sh', 'w')
	f22=open('blast_parser_simple.sh', 'w')
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
			if dedup: fqin = fid+'.dup'
			else: fqin = fid+'.fil'
			f.write(serverTag + ' '+dirscr+'dedup.py '+fqfil+' '+fid+'.dup &\n')
			f1.write(serverTag+' '+dirscr+'fq2faID.py '+fqin+' '+index+' '+fid+'.fa &\n')
			f11.write(serverTag+' '+dirscr+'fq2faID.py '+fastqfile+' '+index+' '+fid+'.fa &\n')
			f2.write(serverTag + ' '+blastnpath+' -task blastn -evalue 1  -max_target_seqs 100000000 -outfmt \'"6  qseqid  sseqid evalue qstart qend sstart send"\'')
			f2.write(' -query '+adaptors+ ' -num_threads '+thread+'  -db '+bpath+index+'_blastdb'+' -out '+fid+'.tab & \n')
			f22.write(serverTag+' '+dirscr+'blast_parser_simple.py '+fid+'.xml '+fid+'.HIV_hits.fa & \n')
			f3.write(serverTag + ' '+dirscr+'blast_trim.py '+fqin+' '+fid+'.tab '+fid+'.ada '+ base+'_'+key+ ' '+str(i)+' &\n') #also recode sequence ID to number
			f4.write(serverTag + ' '+dirscr+'trim_quality.py '+fid+'.ada '+fid+'.trim 33 '+fastqfile+' &\n') #quality trimming
			trimfiles.append(fid+'.trim')
			seqfiles.append(fid+'_sequence.txt')
			if job%nservers==0: 
				f.write('wait\n')
				f1.write('wait\n')
				f11.write('wait\n')
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
	f11.close()
	f22.close()

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
	#clear viralmetagenomics.net #find . -type d -name blast -exec rm -rf "{}" \;
	f2.write('find . -type f ! -name "*.gz"  -delete && find . -type d -empty -delete && rm -rf abyss* soap_out*\n');
	f2.close()

def polyA():
	job=0
	f = open('polyA_raw.sh', 'w')
	f2 = open('polyA_clean.sh', 'w')
	f3= open('plot_poly.sh', 'w')
	for (key, fqfiles) in seeds.items():
		i=1
		for fqfile in fqfiles:
			serverTag = 'ssh '+servers[job%nservers]
			keypath=wd+'/'+base+'/hist/'+key+'_'+str(i)
			
			histR=keypath+'_'+'raw'+'_hist.R'
			outfile=keypath+'_'+'raw'+'_hist.Rout'

			f.write(serverTag+' '+ dirscr+'polyA.py '+wd+'/fastq/'+fqfile+' '+keypath+' raw  &\n')
			f3.write('R CMD BATCH --quiet --vanilla '+histR+' '+outfile+'\n')
			f2.write(serverTag+' '+dirscr+'polyA.py '+wd+'/fastq/'+ key+'_'+str(i)+'_sequence.txt '+keypath+' clean &\n')
			histR=keypath+'_'+'clean'+'_hist.R'
			outfile=keypath+'_'+'clean'+'_hist.Rout'
			f3.write('R CMD BATCH --quiet --vanilla '+histR+' '+outfile+'\n')
			job+=1
			if job%nservers==0: 
				f.write('wait\n')
				f2.write('wait\n')
			i+=1
	f.close()
	f2.close()
	f3.close()



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

def bowtieNT():
	f = open('bowtieNT.sh', 'w')
	f2 = open('bowtiesam2fqNT.sh', 'w')
	f3 = open('bowtieHTML.sh', 'w')
	try: os.mkdir('NT')
	except: pass
	try: os.mkdir(wd+'/'+base)
	except: pass
	try: os.mkdir(wd+'/'+base+'/clark')
	except: pass
	try: os.mkdir(wd+'/'+base+'/clark/fasta')
	except: pass
	
	job=0
	job2=0
	bowindex=bowtieNTIndex
	print bowindex
	csvfiles=[]
	countfiles=[]
	for (key, fqfiles) in seeds.items():
		fafile=wd+'/fastq/'+key+'.fa'
		i=0
		for bowtieind in bowindex:
			job+=1
			i+=1
			sam=wd+'/NT/'+ key+'_'+str(i)+'.sam'
			serverTag = 'ssh '+servers[job%nservers]
			f.write(serverTag+' '+bowtiepath+' --quiet --local --very-fast-local --no-hd --reorder -p 8 -x '+bowtieind+' -f '+fafile+' -S '+sam+' &\n')
			if job%nservers==0: f.write('wait\n')
		job2+=1
		serverTag2 = 'ssh '+servers[job2%nservers]
		start=1
		end=len(bowindex)
		f2.write(serverTag2 + ' '+dirscr+'samNT.py '+key +' '+wd+' '+base+' '+str(start)+' '+str(end)+' & \n')
		f3.write(dirscr+'clark_html.py '+wd+'/'+base+'/clark/'+key+'.count '+wd+'/'+base+'/clark/'+key+'.count.csv '+wd+'/'+base+'/clark/'+key+'.html \n')
		csvfiles.append(wd+'/'+base+'/clark/'+key+'.count.csv')
		countfiles.append(wd+'/'+base+'/clark/'+key+'.count')
		if job2%nservers==0:  f2.write('wait\n')
	f3.write(dirscr+'summaryCount.py '+' '.join(countfiles)+' '+wd+'/'+base+'/clark/all.count\n')
	f3.write(dirscr+'summaryCount.py '+' '.join(csvfiles)+' '+wd+'/'+base+'/clark/all.count.csv\n')
	key='all'
	f3.write(dirscr+'clark_html.py '+wd+'/'+base+'/clark/'+key+'.count '+wd+'/'+base+'/clark/'+key+'.count.csv '+wd+'/'+base+'/clark/'+key+'.html \n')
	f.close()
	f2.close()
	f3.close()

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
		if not keep_human and keep_bac : start, end = '0', '0' #remove human
		elif not keep_human and not keep_bac : start, end = '0', str(len(bowtiebacs)) #remove human
		elif keep_human : start, end = '1', str(len(bowtiebacs))
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
			# # f.write('/mnt/cluster2/xdeng/tools/bwa-0.5.9/bwa aln -n 0.2 -o 3 -t 7 /mnt/cluster/xdeng/hg19/mrnadna.fa fastq/' +index +'_sequence.txt.tmp > fastq/'+ index+'.sai \n')
			# # f.write('/mnt/cluster2/xdeng/tools/bwa-0.5.9/bwa samse /mnt/cluster/xdeng/hg19/mrnadna.fa fastq/' + index+'.sai fastq/'+ index+'_sequence.txt.tmp > fastq/'+ index+'.sam \n')
			# i+=1
		# serverTag2 = 'ssh '+servers[job2%nservers]
		# if pair: f2.write(serverTag2 + ' '+dirscr+'sam2fq.py '+sams[0]+' '+sams[1]+' '+fqfils[0]+' '+fqfils[1]+ ' & \n')
		# else: f2.write(serverTag2 + ' '+dirscr+'sam2fq.py '+sams[0]+' '+fqfils[0]+ ' & \n')
		# job2+=1
		# if job2%nservers==0: f2.write('wait\n')
	# f.close()
	# f2.close()

def bowtieBac(pair, keep_human, keep_bac):
	f = open('bowtieBac.txt', 'w')
	job=0
	if keep_human and not keep_bac: bowindex=bowtiebacs
	elif not keep_human and keep_bac: bowindex= [bowtieindexpath]
	elif not keep_human and not keep_bac: bowindex= [bowtieindexpath]; bowindex.extend(bowtiebacs)
	print '===================================================='
	print bowindex
	print '===================================================='
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
				serverTag=''
				f.write(serverTag+' '+bowtiepath+' --quiet --local --very-fast-local --no-hd --reorder -p 7 -x '+bowtieind+' -U '+fastqfile+' -S '+sam+' \n')
				job+=1
				#if job%nservers==0: f.write('wait\n')
				i+=1
			bj+=1
	#f.write('exit')
	f.close()

def prep_sra():
	kmers = [''.join(p) for p in itertools.product("ACGT", repeat=6)]
	f = open('sra.sh', 'w')
	f2 = open('sra.table', 'w')
	i=0
	for (key, fqfiles) in seeds.items():
		i+=1
		sams=[]
		fqfils=[]
		j=0
		for fqfile in fqfiles:
			j+=1
			fastqfile=wd+'/fastq/'+fqfile
			tag=str(j)+':N:0:'+kmers[i]
			f.write(dirscr+'sra.py '+fastqfile+' '+tag+' '+wd+'/sra.fq.gz\n')
			if j==1: f2.write(fqfile+'\t'+tag)
	f.close()
	f2.close()


def prep_reads(n, length, pair): #filter length fq file and contig, and then combine
	f = open(wd+'/prep_reads.sh', 'w')
	job=0 #rename for Spade
	for (key, vdict) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		f.write(serverTag+' cat '+wd+'/fastq/'+key+'_1_sequence.txt  > '+wd+'/fastq/'+key+'.1.fq &\n')
		if pair: f.write(serverTag+' cat '+wd+'/fastq/'+key+'_2_sequence.txt  > '+wd+'/fastq/'+key+'.2.fq &\n')
		job+=1
		if job%nservers==0: f.write('wait\n')
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
	fk=open('prep_extend.sh', 'w')
	for (key, vdict) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		f.write(serverTag+' '+dirscr+'fqLenFilter.py '+ wd+'/fastq/'+key+'.fq '+' '+wd+'/fastq/'+key+'abyss.fq 35 &\n')
		fk.write(serverTag+' '+dirscr+'fqLenFilter.py '+ wd+'/fastq/'+key+'.fq '+' '+wd+'/fastq/'+key+'.extend.fq 80 &\n')
		job+=1
		if job%nservers==0: f.write('wait\n')
	fk.close()
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
	# job=0
	# for (key, fqfiles) in seeds.items():
		# serverTag = 'ssh '+servers[job%nservers]
		# i=1
		# for fqfile in fqfiles:
			# index=key+'_'+str(i)
			# # raw=wd+'/fastq/'+index+'raw.fq'
			# seqFile = wd+'/fastq/'+index+'_sequence.txt'
			# fqout=wd+'/fastq/'+index+'mira.fastq'
			# f.write(serverTag+' '+dirscr+'mira_shortenID.py '+ seqFile+' '+fqout+' '+str(i)+':N:0:12 &\n')
			# #f.write(serverTag+' '+dirscr+'mira_shortenID.py '+ raw+' '+fqout+' &\n')
			# i+=1
			# job+=1
			# if job%nservers==0: f.write('wait\n')
	# f.close()

def meta_velvet():
	f = open('meta_velvet.sh', 'w')
	# velveth out-dir 51 -fastq -shortPaired HMP.small/SRR041654_shuffled.fastq  HMP.small/SRR041655_shuffled.fastq 
	# velvetg out-dir -exp_cov auto -ins_length 260
	# meta-velvetg out-dir -ins_length 260 | tee logfile
	rval=[]
	job=0
	for (key, vdict) in seeds.items():
		outdir = wd+'/velvet_'+key+'/'
		#try: os.mkdir(outdir)
		#except: pass
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
		ff.write('max_rd_len=500\n[LIB]\n')
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

def abyss(): #filter length fq file and contig, and then combine
	f = open(wd+'/abyss.sh', 'w')
	rval=[]
	job=0
	for (key, vdict) in seeds.items():
		try: os.mkdir(wd+'/abyss_'+key)
		except: pass
		serverTag = '{ time ssh '+servers[job%nservers]
		f.write(serverTag+' "cd '+wd+' && '+abysspath+' -C '+wd+'/abyss_'+key+' name='+key+' k='+abysskmer+' se='+wd+'/fastq/'+key+'abyss.fq" ; } 2> '+wd+'/'+key+'_abyss.time &\n')
		rval.append(wd+'/abyss_'+key+'/'+key+'-unitigs.fa')
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
	rval=[]
	job =0
	for (key, vdict) in seeds.items():
		rv=[]
		for j in xrange(20): #assume 100 chunks, unknown file size but some may not exist
			try: os.mkdir(wd+'/abyss_'+key+'_chunk'+str(j))
			except: pass
			serverTag = 'ssh '+servers[job%nservers]
			f.write(serverTag+' "cd '+wd+' && '+abysspath+' -C '+wd+'/abyss_'+key+'_chunk'+str(j)+' name='+key+'_chunk'+str(j)+ \
			' k='+abysskmer+' se='+wd+'/fastq/'+key+'abyss.fq_'+str(j)+'" &\n')
			#f.write('abyss-pe name='+key+'_chunk'+str(j)+' k='+abysskmer+' se='+wd+'/fastq/'+key+'.fq_'+str(j)+'  \n')
			rv.append(wd+'/abyss_'+key+'_chunk'+str(j)+'/'+key+'_chunk'+str(j)+'-unitigs.fa')
			#rval.append(wd+'/abyss_'+key+'/'+key+'-unitigs.fa')
			job+=1
			if job%nservers==0: f.write('wait\n')
		f1.write('cat '+' '.join(rv) + ' > '+wd+'/'+key+'_abyss_partition.fa\n')
		rval.append(wd+'/'+key+'_abyss_partition.fa')
	f.close()
	f1.close()
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
		f.write(serverTag+' "cd '+wd+' && mira  -t 8 '+wd+'/'+key+'.conf" ; } 2> '+wd+'/'+key+'_mira.time &\n')
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
		if 'T' in assembly_para: combfile.append(contig_out+key+'_T.contig') #trinity
		if 'M' in assembly_para: combfile.append(contig_out+key+'_M.contig') 
		if 's' in assembly_para: combfile.append(contig_out+key+'_s.contig') 
		if 'v' in assembly_para: combfile.append(contig_out+key+'_v.contig') 
		if 'a' in assembly_para: combfile.append(contig_out+key+'_a.contig') 
		f.write('cat '+' '.join(combfile)+ ' > '+wd+'/fastq/'+key+'_contig\n')
		combfile=[]
		f.write(dirscr+'faLenFilter.py '+wd+'/fastq/'+key+'_contig '+wd+'/fastq/'+key+'_contig2 '+str(contigLength1) +'\n')
	f.close()
#move contig for individual assemblers
def moveContig1(assembly_para, r1,r2,r3,r4, r5,r6,r7,r10):
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
		if 'T' in assembly_para:  f.write('mv '+r10[i]+' '+contig_out+key+'_T.contig\n') #velvet
		if 'A' in assembly_para: 
			f.write('cp --dereference '+r3[i]+' '+contig_out+key+'_A.contig\n')
			#f.write('rm '+r3[i]+'\n')#Abyss
		if 'M' in assembly_para: f.write('mv '+r4[i]+' '+contig_out+key+'_M.contig\n') #Mira
		if 's' in assembly_para: f.write('mv '+r5[i]+' '+contig_out+key+'_s.contig\n')
		if 'v' in assembly_para: f.write('mv '+r6[i]+' '+contig_out+key+'_v.contig\n')
		if 'a' in assembly_para: f.write('mv '+r7[i]+' '+contig_out+key+'_a.contig\n')
		# if 'S' in assembly_para: f1.write(dirscr+'statContigs.py '+contig_out+key+'_S.contig  soap '+key+' '+assembly_para+'\n')
		# if 'V' in assembly_para: f1.write(dirscr+'statContigs.py '+contig_out+key+'_V.contig  meta '+key+' '+assembly_para+'\n')
		# if 'A' in assembly_para: f1.write(dirscr+'statContigs.py '+contig_out+key+'_A.contig  abyss '+key+' '+assembly_para+'\n')
		# if 'M' in assembly_para: f1.write(dirscr+'statContigs.py '+contig_out+key+'_M.contig  mira '+key+' '+assembly_para+'\n')
		# if 'S' in assembly_para: f2.write(dirscr+'blastContig.py '+contig_out+key+'_S.contig  soap_'+soapkmer+'_'+key+' '+blastRef+'\n')
		# if 'V' in assembly_para: f2.write(dirscr+'blastContig.py '+contig_out+key+'_V.contig  meta_'+metakmer+'_'+key+' '+blastRef+'\n')
		# if 'A' in assembly_para: f2.write(dirscr+'blastContig.py '+contig_out+key+'_A.contig  abyss_'+abysskmer+'_'+key+' '+blastRef+'\n')
		# if 'M' in assembly_para: f2.write(dirscr+'blastContig.py '+contig_out+key+'_M.contig  mira_'+key+' '+blastRef+'\n')
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
		if 'C' in assembly_para: f.write('mv '+wd+'/fastq/'+key+'_contig3 '+contig_out+key+'_C.contig\n')
		if 'O' in assembly_para: f.write('mv '+wd+'/fastq/'+key+'_contig2-contigs.fa '+contig_out+key+'_O.contig\n')
		# if 'C' in assembly_para: f1.write(dirscr+'statContigs.py '+contig_out+key+'_C.contig cap3 '+key+' '+assembly_para+'\n')
		# if 'O' in assembly_para: f1.write(dirscr+'statContigs.py '+contig_out+key+'_O.contig minimo '+key+' '+assembly_para+'\n')
		# if 'C' in assembly_para: f2.write(dirscr+'blastContig.py '+contig_out+key+'_C.contig cap3_'+key+' '+blastRef+'\n')
		# if 'O' in assembly_para: f2.write(dirscr+'blastContig.py '+contig_out+key+'_O.contig minimo_'+key+' '+blastRef+'\n')
		i+=1
	f1.close()
	f2.close()
	# f.write('mv '+wd+'/assembly.sh '+contig_out+'\n')
	f.write('mv '+wd+'/*.time '+contig_out+'\n') 
	# f.write('rm -rf '+wd+'/fastq/*contig* \nwait\n')
	f.close()



def cap3():
#ssh bsidna3 cap3 /mnt/san2/xdeng/Vhunt/131210_TN_CDC5/fastq/Phan491DNA_contig2 &
	f = open('cap3.sh', 'w')
	job=0
	for (key, vdict) in seeds.items():
		serverTag = '{ time ssh '+servers[job%nservers]
		#f.write(serverTag+' "'+cap3path+' '+wd+'/fastq/'+key+'_contig2" ; } 2> '+wd+'/'+key+'_cap3.time &\n')
		  # -o  N  specify overlap length cutoff > 15 (40)
  # -p  N  specify overlap percent identity cutoff N > 65 (90)
		f.write(serverTag+' "'+cap3path+' '+wd+'/fastq/'+key+'_contig2 -o 25 -p 80" ; } 2> '+wd+'/'+key+'_cap3.time &\n')
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

def minimoFilter(contigLength2):
	f = open('minimoFilter.sh', 'w')
	job=0
	for (key, vdict) in seeds.items():
		#serverTag = 'ssh '+servers[job%nservers]
		serverTag=''
		f.write(serverTag+' '+dirscr+'faLenFilter.py '+wd+'/fastq/'+key+'_contig2-contigs.fa '+wd+'/fastq/'+key+'_contig4 '+str(contigLength2)+' '+base+'_'+key+'\n')
		# job+=1
		# if job%nservers==0: f.write('wait\n')
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
def moveTrinity(contigLength2):
	f = open('moveTrinity.sh', 'w')
	job=0
	for (key, vdict) in seeds.items():
		#serverTag = 'ssh '+servers[job%nservers]
		serverTag=''
		#f.write(serverTag+' cat '+wd+'/fastq/'+key+'_contig2.cap.singlets '+wd+'/fastq/'+key+'_contig2.cap.contigs > '+ wd+'/fastq/'+key+'_contig3\n')
		f.write(serverTag+' '+dirscr+'faLenFilter.py '+wd+'/contig_individual/'+key+'_T.contig '+ wd+'/fastq/'+key+'_contig4 '+str(contigLength2)+' '+base+'_'+key+'\n')
		# job+=1
		# if job%nservers==0: f.write('wait\n')
	f.close()

def reBlastn():
	f = open('blastn1.sh', 'w')
	f2 = open('blastn2.sh', 'w')
	f3 = open('blastn3.sh', 'w')
	f32 = open('blastn32.sh', 'w')
	f33 = open('blastx33.sh', 'w')
	f4 = open('blastn4.sh', 'w')
	f5 = open('blastn_all.sh', 'w')
	try: os.mkdir(wd+'/blastn/')
	except: pass
	ref='vene_ref.fasta'
	os.system('makeblastdb -in '+ref+' -dbtype nucl -parse_seqids -out reference')
	os.system('makeblastdb -in ref_prot.fasta -dbtype prot -parse_seqids -out reference')
	job=0
	outs, outs2, outs3=[], [], []
	for (key, fqfiles) in seeds.items():
		serverTag = 'ssh '+servers[job%nservers]
		fastqfile1=wd+'/fastq/'+fqfiles[0]
		fastqfile2=wd+'/fastq/'+fqfiles[1]
		f.write(serverTag+' '+dirscr+'fq2fa.py '+fastqfile1+' '+fastqfile1+'.fa &\n')
		f.write(serverTag+' '+dirscr+'fq2fa.py '+fastqfile2+' '+fastqfile2+'.fa &\n')
		f2.write(serverTag+' cat '+fastqfile1+'.fa ' +fastqfile2+'.fa > '+wd+'/blastn/'+key+'_reblast.fa &\n')
		out=wd+'/blastn/'+ key+'.out'
		out2=wd+'/blastn/'+ key+'.out2'
		out3=wd+'/blastn/'+ key+'.out3'
		f3.write(serverTag+' '+blastnpath+' -num_threads '+thread+' -outfmt 6 -dbsize 10000000 -searchsp 1000000000 -max_target_seqs 1 -best_hit_overhang 0.1 -best_hit_score_edge 0.1'+' -db '+wd+'/reference'+' -query '+wd+'/blastn/'+key+'_reblast.fa'+' -out '+out+' &\n')
		
		f32.write(serverTag+' '+blastnpath+' -num_threads '+thread+' -perc_identity 100 -ungapped -outfmt 6 -dbsize 10000000 -searchsp 1000000000 -max_target_seqs 1 -best_hit_overhang 0.1 -best_hit_score_edge 0.1'+' -db '+wd+'/reference'+' -query '+wd+'/blastn/'+key+'_reblast.fa'+' -out '+out2+' &\n')
		f33.write(serverTag+' '+blastxpath+' -num_threads '+thread+' -outfmt 6 -dbsize 10000000 -searchsp 1000000000 -max_target_seqs 1 -best_hit_overhang 0.1 -best_hit_score_edge 0.1'+' -db '+wd+'/reference'+' -query '+wd+'/blastn/'+key+'_reblast.fa'+' -out '+out3+' &\n')
		outs.append(out)
		outs2.append(out2)
		outs3.append(out3)
		job+=1
		if job%nservers==0: f.write('wait\n')
	f.close()
	f2.close()
	f3.close()
	f32.close()
	f33.close()
	#f4.write(dirscr+'tally_blastn.py '+ref+' '+' '.join(outs)+' '+wd+'/blastn_tally.txt\n')
	f4.write(dirscr+'tally_blastn.py '+ref+' '+' '.join(outs2)+' '+wd+'/blastn_tally_perfectmatch.txt\n')
	#f4.write(dirscr+'tally_blastn.py '+ref+' '+' '.join(outs3)+' '+wd+'/blastx_tally_prot.txt\n')
	f4.close()
	f5.write('source blastn1.sh\nwait\n')
	f5.write('source blastn2.sh\nwait\n')
	#f5.write('source blastn3.sh\nwait\n')
	f5.write('source blastn32.sh\nwait\n')
	#f5.write('source blastx33.sh\nwait\n')
	f5.write('source blastn4.sh\n')
	f5.close()

def reAssemble(fafile):
	f0 = open('reAssembleAll.sh', 'w')
	f0.write('source reAssemble.sh\nwait\n')
	f0.write('source sam2count.sh\nwait\n')
	f0.write('source annotate_reAssemble.sh\nwait\n')
	f0.close()
	f = open('reAssemble.sh', 'w')
	try: os.mkdir(wd+'/bowtieContig/')
	except: pass
	try: os.mkdir(wd+'/'+base+'/')
	except: pass
	try: os.mkdir(wd+'/'+base+'/reAssemble/')
	except: pass
	if fafile=='': 
		f.write('cat ')
		for (key, vdict) in seeds.items():
			f.write( wd+'/fastq/'+key+'_contig4 ')
		f.write(' > '+wd+'/bowtieContig/all.cap3\n')
		f.write(cap3path+' '+wd+'/bowtieContig/all.cap3 -o 25 -p 85 -h 40\n')
		f.write(' cat '+wd+'/bowtieContig/all.cap3.cap.singlets '+wd+'/bowtieContig/all.cap3.cap.contigs > '+ wd+'/bowtieContig/all.cap3cap3\n')
		fafile=wd+'/fastq/reAssemble_c'
		f.write(dirscr+'faLenFilter.py '+wd+'/bowtieContig/all.cap3cap3 '+ fafile+' 1000 dummy\n')
	# else: fafile=wd+'/beatrixContig.fa'
	f.write('/mnt/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build '+fafile+' '+wd+'/bowtieContig/allbow\n')

	job=0
	bowindex= wd+'/bowtieContig/allbow'
	for (key, fqfiles) in seeds.items():
		fastqfile1=wd+'/fastq/'+fqfiles[0]
		fastqfile2=wd+'/fastq/'+fqfiles[1]
		sam=wd+'/bowtieContig/'+ key+'.sam'
		serverTag = 'ssh '+servers[job%nservers]
		# at least 30 bp perfect match -L 30, default is 20 
		f.write(serverTag+' '+bowtiepath+' --fast-local --quiet --local --no-hd --reorder -p 7 -L 30 -x '+bowindex+' -1 '+fastqfile1+' -2 '+fastqfile2+' -S '+sam+' &\n')
		job+=1
		if job%nservers==0: f.write('wait\n')
	f2 = open('sam2count.sh', 'w')
	job=0
	countfiles=[]
	for (key, fqfiles) in seeds.items():
		sam=wd+'/bowtieContig/'+ key+'.sam'
		count=wd+'/bowtieContig/'+ key+'.count'
		serverTag = 'ssh '+servers[job%nservers]
		f2.write(serverTag+' '+dirscr+'sam2count.py '+sam+' '+count+' &\n')
		countfiles.append(count)
		job+=1
		if job%nservers==0: f2.write('wait\n')
	f2.write('wait\n')
	f2.write(dirscr+'tally.py '+fafile+' '+' '.join(countfiles)+' '+wd+'/'+base+'/reAssemble/tally.txt\n')
	f.close()
	f2.close()
	
	f3=open('annotate_reAssemble.sh', 'w')
	f3.write('cat '+wd+'/'+base+'/fasta/reAssemble_*.fa > '+wd+'/'+base+'/reAssemble/hits0.fa \n')
	f3.write('cat '+wd+'/'+base+'/mystery/reAssemble_m_complex_sort.fasta > '+wd+'/'+base+'/reAssemble/mystery0.fa \n')
	f3.write(dirscr+'annotate_contig.py '+wd+'/'+base+'/reAssemble/tally.txt '+wd+'/'+base+'/reAssemble/hits0.fa '+ \
			wd+'/'+base+'/reAssemble/mystery0.fa ' + wd+'/'+base+'/reAssemble/reAssemble.txt ' \
			+wd+'/'+base+'/reAssemble/hits.fa '+ wd+'/'+base+'/reAssemble/mystery.fa \n')
	f3.close()

def Assembly_individual(assembly_para):
	partition()
	sf = open('assembly_individual.sh', 'w')
	r1,r2,r3,r4, r5,r6,r7=None, None, None, None, None, None, None
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
	if 'T' in assembly_para: 
		r10=trinity(RAM)
		sf.write('source trinity.sh\nwait\n')
	else: r10=None
	moveContig1(assembly_para, r1,r2,r3,r4, r5,r6,r7, r10)
	sf.write('echo '+password+' |find . "*" -print0 | sudo xargs -0 chmod 777\n')
	sf.write('source moveContig.sh\n')
	contig_out = wd+'/contig_individual/'
	#sf.write('source statContig.sh > '+contig_out+'contig.log\n')
	#sf.write('source blastContig.sh > '+contig_out+'contig.blast\n')
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
		minimoFilter(contigLength2)
		os.system('mv minimo.sh '+contig_out+'\nwait\n')
		os.system('mv minimoFilter.sh '+contig_out+'\nwait\n')
		sf.write('source '+contig_out+'minimo.sh\nwait\n')
		sf.write('source '+contig_out+'minimoFilter.sh\nwait\n')

	moveContig2(assembly_para)
	contig_out = wd+'/contig_'+assembly_para+'/'
	try: os.mkdir(contig_out)
	except: pass
	sf.write('echo '+password+' |find . "*" -print0 | sudo xargs -0 chmod 777\n')
	sf.write('source '+contig_out+'moveContig.sh\n')
	sf.write('source '+contig_out+'statContig.sh > '+contig_out+'contig.log\n')
	sf.write('source '+contig_out+'blastContig.sh > '+contig_out+'contig.blast\n')
	sf.close()

def AddPipe(para, sf):
	assembly_para=para
	Assembly_combine(assembly_para)
	sf.write('source '+wd+'/contig_'+assembly_para+'/'+'assembly_combine.sh\n')

	
#.contig3 has been moved to contig_out directory from fastq directory, .contig4 still there and will be merged to reads
def combineContig_reads(n, skipread): #filter length fq file and contig, and then combine
	f = open(wd+'/combineContig_reads.sh', 'w')
	if doreAssemb: 
		f.write(dirscr+'splitQuery.py '+wd+'/fastq/reAssemble_c '+' '+str(n)+'\n')
	else:
		for (key, vdict) in seeds.items():
			contig = wd+'/fastq/'+key+'_contig4 ' #contig combined with reads
			combine = wd+'/fastq/'+key+'_c ' #contig combined with reads
			if skipread: f.write('cat '+contig+' > '+combine +'\n')
			else: f.write('cat '+contig+' '+wd+'/fastq/'+key+'.fa > '+combine +'\n')
			f.write(dirscr+'splitQuery.py '+combine+' '+str(n)+'\n')
	f.close()
	
def justBlast(n): #filter length fq file and contig, and then combine
	f = open(wd+'/justBlast.sh', 'w')
	for (key, fqfiles) in seeds.items():
		combine = wd+'/fastq/'+key+'_c ' #contig combined with reads
		f.write('cat '+wd+'/fastq/'+fqfiles[0].replace('fastq', 'fasta')+' > '+combine +'\n')
		f.write(dirscr+'splitQuery.py '+combine+' '+str(n)+'\n')
	f.close()

def trinity(RAM): #This is actually spade
	f = open('trinity.sh', 'w')
	rval=[]
	
	job=0
	for (key, fqfiles) in seeds.items():
		fq = [wd+'/fastq/'+fqfile for fqfile in fqfiles]
		outdir = wd+'/trinity_'+key+'/'
		serverTag = '{ time ssh '+servers[job%nservers]
		file1 = wd+'/fastq/'+key+'.1.fq' #fq[0]#
		#file2 = fq[1]#wd+'/fastq/'+key+'_2_sequence.txt'
		if pair:
			file2 = wd+'/fastq/'+key+'.2.fq'	#fq[1]#
			f.write(serverTag+ ' "'+spadepath+'spades.py -m '+ SI[servers[job%nservers]][1]+' --meta -k 21,33,55,77 -t 48 -1 '+file1+' -2 '+file2+' -o '+outdir+'"  ; } 2> '+wd+'/'+key+'_trinity.time &\n')
		else: f.write(serverTag+ ' "'+spadepath+'spades.py -m '+ SI[servers[job%nservers]][1]+' -k 21,33,55,77 --careful -t 48 -s '+file1+' -o '+outdir+'"  ; } 2> '+wd+'/'+key+'_trinity.time &\n')
		# if pair: 
			# #Trinity --seqType fq --JM 100G --left reads_1.fq  --right reads_2.fq --CPU 6
			# f.write(serverTag+' "'+trinitypath+'/Trinity --seqType fq --max_memory '+RAM+' --no_bowtie --output '+outdir+' --left '+file1 + \
					# ' --right '+file2+' --CPU 8 "  ; } 2> '+wd+'/'+key+'_trinity.time &\n')
		# else:
			# f.write(serverTag+' "'+trinitypath+'/Trinity --seqType fq --max_memory '+RAM+' --no_bowtie --output '+outdir+' --single '+file1 + \
					# ' --CPU 8 "  ; } 2> '+wd+'/'+key+'_trinity.time &\n')
		trinityOut=outdir+'/contigs.fasta'
		rval.append(trinityOut)
		job+=1
		if job%nservers==0: f.write('wait\n')
	f.close()
	return rval

def blastVirus(n, hsp):
	cwd = os.path.basename(os.getcwd())
	tt=open('vfam.log', 'w')
	tt.close()
	tt=open('blast_virus.log', 'w')
	tt.close()
	f = open(wd+'/blast_virus.txt', 'w')
	fff=open(wd+'/fixDiamondXML.sh', 'w')
	f1 = open(wd+'/blast_nr.txt', 'w')
	f8 = open(wd+'/rap_nr.txt', 'w')
	f2 = open(wd+'/blast_virus_parser.sh', 'w')
	f3 = open(wd+'/blast_nr_filter.sh', 'w')
	f35 = open(wd+'/diamond_nr_filter.sh', 'w')
	f4 = open(wd+'/blast_output_merge.sh', 'w')
	f40 = open(wd+'/blast_output_sort.sh', 'w')
	f400 = open(wd+'/plot_pie.sh', 'w')
	f5 = open(wd+'/blast_nr_mystery.txt', 'w')
	f6 = open(wd+'/blast_nr_mystery_filter.sh', 'w')
	f7 = open(wd+'/movetowww.sh', 'w')
	f9 = open(wd+'/vfam.sh', 'w')
	f10 = open(wd+'/vfam_annot.sh', 'w')
	f11=open(wd+'/diamond_nr.txt', 'w')
	#f12=open(wd+'/diamond_nr_view.sh', 'w')
	f13=open(wd+'/mergeSig.sh', 'w')
	f50=open(wd+'/dna2prot.sh', 'w')
	f65=open(wd+'/mergeTable.sh', 'w')
	hittable=wd+'/'+base+'/hitTable'
	try: os.mkdir(wd+'/blast_virus_out/')
	except: pass
	try: os.mkdir(wd+'/diamond_nr_out/')
	except: pass
	try: os.mkdir(wd+'/blast_nr_out/')
	except: pass
	try: os.mkdir(wd+'/blast_filter_out/')
	except: pass
	try: os.mkdir(wd+'/'+base+'/')
	except: pass
	try: os.mkdir(wd+'/'+base+'/mystery/')
	except: pass
	try: os.mkdir(wd+'/'+base+'/fastq/')
	except: pass
	try: os.mkdir(wd+'/'+base+'/hist/')
	except: pass
	try: os.mkdir(wd+'/'+base+'/pie/')
	except: pass
	try: os.mkdir(wd+'/blast_filter_out/tmp/')
	except: pass
	# f5 = open(wd+'/blast_filter_out/index.html', 'w')
	# f5.write('<html>')
	dbname= os.path.abspath(virusdbpath)
	dbnameDNA= os.path.abspath(virusdbpathDNA)
	
	dbnr= os.path.abspath(nvnrdbpath)
	job, job2=0,0
	allout = 'cat '
	allmys = 'cat '
	allmysout = 'cat '
	allmysout2 = 'cat '
	allcombine='cat '
	jj=0
	outtable=' '+wd+'/blast_filter_out/table/all_blast_filter.txt '
	for (key, vdict) in seeds.items():
		jj+=1
		queryname = wd+'/fastq/'+key+'_c'
		sigfa=wd+'/fastq/'+key+'_sig'
		# diamondout1=wd+'/diamond_nr_out/'+key
		diamondout2=wd+'/diamond_nr_out/'+key+'.m8'
		mysfile=wd+'/'+base+'/mystery/'+key+'_m.fasta'
		mysoutfile=wd+'/'+base+'/mystery/'+key+'_m.fasta.out'
		mysoutfile2=wd+'/'+base+'/mystery/'+key+'_m.fasta.out2'
		mysfile2=wd+'/'+base+'/mystery/'+key+'_m_complex.fasta'
		mysfile3=wd+'/'+base+'/mystery/'+key+'_m_complex_sort.fasta'
		outname = wd+'/blast_virus_out/'+key+'_blast.xml'
		outname1 = wd+'/blast_nr_out/'+key+'_blast.xml'
		outtxt = wd+'/blast_filter_out/'+key+'_blast_filter.txt'
		outtable+=' '+wd+'/blast_filter_out/table/'+key+'_blast_filter.txt'
		# htmlFile =os.path.basename(os.path.splitext(outtxt)[0]+'.html')
		# f5.write('<a href="'+htmlFile+'">'+htmlFile+'</a><br>')
		merge='cat '
		mys='cat '
		mysout='cat '
		mergeSig='cat '
		serverTag2 = 'ssh '+servers[jj%nservers]
		for i in xrange(n):
			serverTag = 'ssh '+servers[job%nservers]

			subquery=queryname+'_'+str(i)
			sigfaname = subquery+'_s' #output significant virus fa file for NR
			submyspre = wd+'/blast_filter_out/tmp/'+key+'_prem.fasta'+'_'+str(i) #output mysterious contigs
			submys = wd+'/blast_filter_out/tmp/'+key+'_m.fasta'+'_'+str(i) #output mysterious contigs filtered by NR
			submys_p = wd+'/blast_filter_out/tmp/'+key+'_m.fasta'+'_'+str(i)+'_p' #output mysterious contigs filtered by NR
			submysout = wd+'/blast_filter_out/tmp/'+key+'_m.fasta'+'_out'+str(i) #output mysterious contigs filtered by NR
			submysout2 =wd+'/blast_filter_out/tmp/'+key+'_m2.fasta'+'_out'+str(i)
			submysxml = wd+'/blast_filter_out/tmp/'+key+'_m.xml_'+str(i) #output mysterious contigs
			subxml0=outname+'_'+str(i)+'.pre' #these are diamond XML output with bug need fixing
			subxml=outname+'_'+str(i)
			subxml1=outname1+'_'+str(i)
			subout=outtxt+'_'+str(i)

			

			
			if doDNA: 
				f.write(blastnpath+' -num_threads '+thread+'  -evalue '+EVALUE+' -outfmt 5 -dbsize 10000000 -searchsp 1000000000 -max_target_seqs 1')
				f.write(' -db_soft_mask 11 -best_hit_overhang 0.1 -best_hit_score_edge 0.1') #masking and best hit
				f.write(' -db '+dbnameDNA+ ' -query '+ subquery+' -out '+subxml+'\n')
			elif doDiamondOnly:
				f.write(newdiamondpath+' blastx -d '+diamondvirusdb+' --quiet --sensitive --max-target-seqs 1 --outfmt 5 -q '+subquery+' -o '+subxml0+'\n')
				fff.write(serverTag+' '+dirscr+'fixDiamondXML.py '+subxml0+' '+subxml+' & \n')
			else:
				f.write(blastxpath+' -num_threads '+thread+'  -evalue '+EVALUE+' -outfmt 5 -dbsize 10000000 -searchsp 1000000000 -max_target_seqs 1')
				f.write(' -db_soft_mask 21 -best_hit_overhang 0.1 -best_hit_score_edge 0.1') #masking and best hit
				f.write(' -db '+dbname+ ' -query '+ subquery+' -out '+subxml+'\n')
				
			f1.write(blastxpath+' -num_threads '+thread+'  -evalue '+EVALUE+' -outfmt 5 -dbsize 10000000 -searchsp 1000000000 -max_target_seqs 1')
			f1.write(' -best_hit_overhang 0.1 -best_hit_score_edge 0.1') #best hit only
			f1.write(' -db '+dbnr+ ' -query '+ sigfaname+' -out '+subxml1+'\n')
			f8.write(rapsearch+' -q '+sigfaname+' -d '+rapdbnr+' -o '+subxml1+' -z 8 -b 1 -v 1 -x t\n')

			f2.write(serverTag+' '+dirscr+'blast_parser.py '+subxml+ ' '+ subquery+ ' '+ sigfaname+ ' '+submyspre+' '+str(myslen)+' '+EVALUE+' '+wd+'/blast_virus.log & \n')
			f5.write(blastxpath+' -num_threads '+thread+'  -evalue 0.001 -outfmt 5 -dbsize 10000000 -searchsp 1000000000 -max_target_seqs 1')
			f5.write(' -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -db '+dbnr+ ' -query '+ submyspre+' -out '+submysxml+'\n')
			#diamond blastx -d nr -q reads.fna -a matches 
			#diamond view -a matches.daa -o matches.m8
			# f11.write(diamondpath+' blastx --max-target-seqs 1 -d '+diamonddb+ ' -q '+ sigfaname+' -a '+subxml1+'\n')
			# f12.write(serverTag+' '+diamondpath+' view -a '+subxml1+'.daa'+ ' -o '+ subxml1+' &\n')
			f6.write(serverTag+' '+dirscr+'blast_filter_NR_mys.py '+ submysxml +' '+submyspre+' '+ submys+' &\n')
			f50.write(serverTag+' '+dirscr+'dna2prot.py '+ submys +' '+submys_p+' &\n')
			f9.write(serverTag+' '+hmmerpath+' --cpu '+thread+' -E 0.001 -A '+submysout+' --tblout '+submysout2+' '+vfampath+' '+submys_p+ ' > 1 &\n')
			# f10.write(serverTag+' '+dirscr+'hmmer_filter_NR_mys.py '+ submysxml +' '+submyspre+' '+ submys+' &\n')
			if dorapsearch: subxml1+='.xml'
			f3.write(serverTag+' '+dirscr+'blast_filter_NR.py '+ subxml + '  '+ subxml1+'  '+ subquery+' '+subout+' '+EVALUE+' '+hsp+' & \n')

			f35.write(serverTag+' '+dirscr+'diamond_filter_NR.py '+ subxml + '  '+ diamondout2+'  '+ subquery+' '+subout+' '+EVALUE+' '+hsp+' & \n')
			job+=1
			if job%nservers==0: 
				f9.write('wait\n')
				fff.write('wait\n')
				f2.write('wait\n')
			#if job%(nservers*4)==0: #send 4 jobs at a batch
				f3.write('wait\n')
				f35.write('wait\n')
				f6.write('wait\n')
				f50.write('wait\n')
			merge+=subout+' '
			mys+=submys+' '
			mysout+=submysout2+' '+submysout+' '
			mergeSig+=sigfaname+' '
		
		merge+='> '+outtxt
		mys+='> '+mysfile
		mysout+='> '+mysoutfile
		mergeSig+= '> '+sigfa
		f13.write(mergeSig+'\n')
		f11.write(newdiamondpath+' blastx --quiet --max-target-seqs 1 --outfmt 6 -d '+diamonddb+ ' -q '+ sigfa+' -o '+diamondout2+'\n') 
		#f12.write(serverTag2+' '+diamondpath+' view -a '+diamondout1+'.daa'+ ' -o '+ diamondout2+' &\n')
		#f12.write(diamondpath+' view -a '+diamondout1+'.daa'+ ' -o '+ diamondout2+' \n')
		f10.write(serverTag2+' '+dirscr+'hmmer_annot.py '+hmmannot+' '+mysoutfile+' '+mysoutfile2+' '+wd+'/vfam.log &\n')
		allout+= outtxt + ' '
		allmys+= mysfile + ' '
		allmysout+= mysoutfile + ' '
		allmysout2+= mysoutfile2 + ' '
		f4.write(merge+'\n')
		f4.write(mys+'\n')
		f4.write(mysout+'\n')
		combine = wd+'/fastq/'+key+'.fa' 
		if 'reAssemb' not in key: allcombine+= combine +' '
		f40.write(serverTag2+ ' " '+ dirscr+'blast_output_sort.py '+ outtxt +' '+combine+' '+dirscr+' '+base+' '+cwd+' && ')
		f40.write(dustmaskerpath+' -in '+mysfile+' -infmt fasta -outfmt fasta -out '+mysfile2+' && ')
		f40.write(dirscr+'faSort.py '+mysfile2+' '+mysfile3+' True " & \n')
		f400.write(dirscr+'plot_pie.py '+ wd+'/'+base+'/pie/'+key+'_blast_filter.txt\n')
		if jj%nservers==0:
			f10.write('wait\n')
			f40.write('wait\n')
			#f12.write('wait\n')

	f65.write(dirscr+'mergeTable.py '+outtable+' '+hittable+'\n')
	outtxt = wd+'/blast_filter_out/'+'all'+'_blast_filter.txt'
	mysfile =wd+'/'+base+'/mystery/'+'all'+'_m.fasta'
	mysfile =wd+'/'+base+'/mystery/'+'all'+'_m.fasta.out'
	mysfile2=wd+'/'+base+'/mystery/all_m_complex.fasta'
	mysfile3=wd+'/'+base+'/mystery/all_m_complex_sort.fasta'
	mysoutfile=wd+'/'+base+'/mystery/all_m.fasta.out'
	mysoutfile2=wd+'/'+base+'/mystery/all_m.fasta.out2'
	
	
	combine =wd+'/fastq/'+'all'+'.fa'
	allout +='> '+ outtxt
	allmys +='> '+ mysfile
	#allmysout +='> '+ mysoutfile
	allmysout2 +='> '+ mysoutfile2
	allcombine += '> '+ combine
	f4.write(allout+'\n')
	f4.write(allmys+'\n')
	#f4.write(allmysout+'\n')
	f10.write(allmysout2+'\n')
	#f10.write(dirscr+'hmmer_annot.py '+hmmannot+' '+mysoutfile+' '+mysoutfile2+' &\n')
	f4.write(allcombine+'\n')
	#f40.write(dirscr+'blast_output_sort.py '+ outtxt +' '+combine+' '+dirscr+'\n')
	f40.write('wait\n')
	f40.write(serverTag2+' " '+dirscr+'blast_output_sort.py '+ outtxt +' '+'--'   +' '+dirscr+' '+base+' '+cwd+' && ') #skip getting reverse strand
	f40.write(dustmaskerpath+' -in '+mysfile+' -infmt fasta -outfmt fasta -out '+mysfile2+' && ')
	f40.write(dirscr+'faSort.py '+mysfile2+' '+mysfile3+' " &\n')
	f400.write(dirscr+'plot_pie.py '+ wd+'/'+base+'/pie/all_blast_filter.txt\n')
	# f7.write('rm -rf '+wd+'/'+base+'/aln/\n')
	# f7.write('rm -rf '+wd+'/'+base+'/fasta/\n')
	# f7.write('rm -rf '+wd+'/'+base+'/tmp/\n')
	# f7.write('rm -rf '+wd+'/'+base+'/pie/\n')
	f7.write('cp -r '+wd+'/blast_filter_out/aln '+wd+'/'+base+'/\n')
	f7.write('cp -r '+wd+'/blast_filter_out/fasta '+wd+'/'+base+'/\n')
	f7.write('cp -r '+wd+'/blast_filter_out/pie '+wd+'/'+base+'/\n')
	f7.write('cp -r '+wd+'/blast_filter_out/*.html '+wd+'/'+base+'/\n')
	f7.write('cp -r '+wd+'/blast_filter_out/*.xls '+wd+'/'+base+'/\n')
	f7.write('cp -r '+wd+'/*.html '+wd+'/'+base+'/\n')
	f7.write('cp '+dirscr+'*.php '+wd+'/'+base+'/\n')
	if genBowtie: f7.write('cp -r '+wd+'/fastq/*.fastq* '+wd+'/'+base+'/fastq/\n')

	
	f.close()
	f1.close()
	f2.close()
	f3.close()
	f4.close()
	f5.close()
	f6.close()
	f7.close()
	f8.close()
	f9.close()
	f50.close()
	f10.close()
	f40.close()
	f400.close()
	f11.close()
	#f12.close()
	f13.close()
	f65.close()
	f35.close()
	fff.close()



if __name__ == "__main__":
	wd = os.path.abspath(os.path.dirname('.')).replace('san2', 'cluster2')
	stats=defaultdict()
	tcga=False
	if tcga:
		bam2sam()
		#tar2fastq()
		of=open('tcga_prep.sh', 'w')
		of.write('source bam2sam.sh\nwait\n')
		of.write('source samfilterfq.sh\nwait\n')
		#of.write('source tar2fq.sh\nwait\n')
		sys.exit()
	
	password='Welcome39'
	n=50 #50 #number fo splits
	length=50 #lengh reads for blast
	contigLength1 = 300 #length before cap3
	contigLength2 = 1500 #length after cap3 before blast
	myslen = 1000
	fasta=False
	bam=False
	mira_local=False
	sra=False
	keep_human=False
	keep_bac=True
	pair=False
	mergePair=False
	phage='False' #'Both' #'False', 'True'
	rm_adaptor=False
	skipread=False
	doMyth=False
	doAssembly= 'no' # 'denovo'#, 'trinity' #, 'no', #Trinity is spade
	#doAssembly=  #'trinity' #, 'no', #Trinity is spade
	dedup=True
	doHmmer=False
	doreAssemb=False
	internal=False
	dorapsearch=False
	donrfilter=True
	doNT=False
	genBlast=True
	genBowtie=True 
	doClark=False
	virusdbpath='/mnt/cluster/xdeng/blastdb/virus_mask'
	doDiamondOnly=False
	doDNA=False
	if doDiamondOnly:
		n=1
	#virusdbpath='/mnt/cluster/xdeng/blastdb/viral_pasteur_mask'

	thread = '48'
	EVALUE='0.01'#'0.01' #threshold for virus
	abysskmer='31'
	soapkmer='31'
	metakmer='31'
	hsp='NO' # show whole reads, if 'YES' show only blast hsp
	
	#genseedfile()
	# if pair:
	#stats,seeds=mergeAB(reAssemble='no')
	
	stats,seeds=readSeeds2()
	# else:
		# seeds=readSeeds1()
	if not wd.startswith('/mnt/'): wd='/mnt'+wd
	base = os.path.basename(wd)
	print 'path',wd
	print 'thread', thread
	print 'base', base
	print 'numBarcodes', len(seeds)
	for key in seeds.keys():
		print len(seeds[key]), '--', key, seeds[key]
	print 'bam', bam
	print 'dedup', dedup
	print 'fasta', fasta
	print 'keep_human', keep_human
	print 'keep_bac', keep_bac
	print 'pairend', pair
	print 'rm_adaptor', rm_adaptor
	print 'skipread', skipread
	print 'hsp', hsp
	print 'metakmer', metakmer
	print 'soapkmer', soapkmer
	print 'abysskmer', abysskmer
	print 'servers', servers
	print 'EVALUE', EVALUE
	print 'doMyth', doMyth
	print 'doAssembly', doAssembly
	print 'doHmmer', doHmmer
	print 'read threshold', length
	print 'contig before cap3', contigLength1
	print 'contig after cap3', contigLength2
	print 'mys threshold', myslen
	print 'internal', internal
	print 'phage', phage
	print 'rapsearch', dorapsearch
	print 'nrfilter', donrfilter
	print 'n=', n
	print 'length', length
	print 'bowtiebacs', bowtiebacs
	print 'doreAssemb', doreAssemb
	print 'password', password
	print 'genBlast', genBlast
	print 'genBowtie', genBowtie
	print 'doDiamondOnly', doDiamondOnly
	print 'pair=', pair
	print 'clark=',doClark
	print 'doDNA=', doDNA
	if phage=='True': #virusdbpath=phan_phagedbpath 
		virusdbpath=phagedbpath
	elif phage=='False': #virusdbpath=phan_phagedbpath 
		virusdbpath=virusdbpath
	elif phage=='Both':
		virusdbpath=virusphagedbpath
	print 'virusdbpath', virusdbpath

	#zip2fq()
	# sys.exit()
	#renamefastq()
	fff=open(wd+'/server.txt', 'w')
	for s in servers:
		fff.write(s+'\n')
	fff.close()
	oof=open(wd+'/run.log', 'w')
	for key in seeds.keys():
		print >> oof, key, seeds[key]
	print >>oof, 'bam', bam
	print >>oof, 'dedup', dedup
	print >>oof, 'keep_human', keep_human
	print >>oof, 'keep_bac', keep_bac
	print >>oof,'pairend', pair
	print >>oof,'rm_adaptor', rm_adaptor
	print >>oof,'skipread', skipread
	print >>oof,'hsp', hsp
	print >>oof, 'metakmer', metakmer
	print >>oof, 'soapkmer', soapkmer
	print >>oof, 'abysskmer', abysskmer
	print >>oof, 'servers', servers
	print >>oof, 'EVALUE', EVALUE
	print >>oof, 'doMyth', doMyth
	print >>oof, 'doAssembly', doAssembly
	print >>oof, 'read threshold', length
	print >>oof, 'contig before cap3', contigLength1
	print >>oof, 'contig after cap3', contigLength2
	print >>oof, 'mys threshold', myslen
	print >>oof, 'virusdbpath', virusdbpath
	print >>oof, 'numSplit=', n
	oof.close()
	# sampleFastq()
	# sys.exit()

	if mergePair: mergePairFq()
	# if pair: 
		# seeds, pair = mergePairFq()
	# print seeds, pair

	prep_sra()
	polyA()
	if fasta: fa2fq()
	if bam: bam2fq()
	# bowtieHuman(pair)
	bowtieNT()
	skipbowtie()
	bowtieBac(pair, keep_human, keep_bac)
	sam2fq(pair, keep_human)
	justBlast(n)
	check_pair_overlap()
	prepBlastFile_adaptor()
	trim()
	trinity ('')
	skip_adaptor()
	fq_check()
	prep_reads(n, length, pair)
	if doreAssemb: 
		#reAssemble()
		reAssemble( wd+'/ref.fasta')
		seeds={}
		seeds['reAssemble']=[]
	combineContig_reads(n, skipread)
	blastVirus(n, hsp)
	#hmmer(n,hsp)
	clean_dir()
	prepBlastFile()
	Clark()
	#reBlastn()

	#if pair: prepPriceFile() #not doing price

	sf=open('pipeline_run.sh', 'w')
	
	if fasta: sf.write('source fa2fq.sh\nwait\n')
	if bam: sf.write('source bam2fq.sh\nwait\n')
	# os.system('PriceTI  -fp s1.fq s2.fq 300 -icf seeds.fa 1 1 5  -nc 30 -dbmax 72 -mol 35 -tol 20 -mpi 80 -target 90 2 2 2 -a 7 -o price.fa')
	if mergePair: sf.write('source mergePairFq.sh\nwait\n')
	if sra: 
		sf.write('source sra.sh\n')
	if keep_human and keep_bac:
		sf.write('source skipbowtie.sh\nwait\n')
	else:
		sf.write(dirscr+'schedule2.py bowtieBac.txt server.txt\nwait\n')
		#sf.write('source bowtieBac.sh\nwait\n')
		sf.write('source bowtiesam2fq.sh  >bowtiesam2fq.log\nwait\n')
	if dedup:
		sf.write('source clonetrim.sh >clonetrim.log \nwait\n')
	if rm_adaptor:
		sf.write('source clonefq2fa.sh \nwait\n')
		sf.write('source prepBlastFile_adaptor.sh \n')
		sf.write('source blast_adaptor.sh \nwait\n')
		sf.write('source blasttrim.sh >blasttrim.log \nwait\n')
		sf.write('source qualitytrim.sh >qualitytrim.log \nwait\n')
	else: sf.write('source skipadaptor.sh \nwait\n')
	sf.write('source fq_clean.sh \nwait\n')
	#sf.write('source check_pair.sh >stats.log \nwait\n')
	sf.write('source polyA_raw.sh >raw.log \nwait\n')
	sf.write('source polyA_clean.sh >clean.log \nwait\n')
	sf.write('echo '+password+' |find . "*" -print0 | sudo xargs -0 chmod 777\n')
	
	sf.write('source prep_reads.sh\nwait\n')
	if doClark:
		sf.write('source clark.sh\nwait\n')
		sf.write('source clark_process.sh\nwait\n')
	if doNT:
		sf.write('source bowtieNT.sh\nwait\n')
		sf.write('source bowtiesam2fqNT.sh\nwait\n')
		sf.write('source bowtieHTML.sh\nwait\n')
	
	if doAssembly=='denovo':
		assembly_para='SAVa' #SAOP, Abyss, Velvet, Mira, Cap3, no Minimo
		Assembly_individual(assembly_para)
		sf.write('source assembly_individual.sh\n')
		AddPipe('SAVaC', sf)
	elif doAssembly=='trinity':
		assembly_para='T' #SAOP, Abyss, Velvet, Mira, Cap3, no Minimo
		Assembly_individual(assembly_para)
		sf.write('source assembly_individual.sh\nwait\n')
		moveTrinity(contigLength2)
		sf.write('source moveTrinity.sh\nwait\n')

	sf.write('source combineContig_reads.sh >combine.log \nwait\n')
	#sf.write('source justBlast.sh\n')
	if doHmmer:
		sf.write('source dna2prot.sh \nwait\n')
		sf.write(dirscr+'schedule2.py hmmer_virus.txt server.txt 4\nwait\n')
		sf.write('source hmmer_virus_parser.sh >hmmer_virus.log  \nwait\n')
	else: 
		sf.write(dirscr+'schedule2.py blast_virus.txt server.txt\nwait\n')
		if doDiamondOnly:
			sf.write('source fixDiamondXML.sh\nwait\n')
		sf.write('source blast_virus_parser.sh \nwait\n')
	if doMyth:
		if doHmmer:
			sf.write('source dna2protmys.sh \nwait\n')
			sf.write(dirscr+'schedule2.py hmmer_nr_mystery.txt server.txt 4\nwait\n')
			sf.write('source hmmer_nr_mystery_filter.sh \nwait\n')
		else:
			sf.write(dirscr+'schedule2.py blast_nr_mystery.txt server.txt\nwait\n')
			sf.write('source blast_nr_mystery_filter.sh \nwait\n')
			sf.write('source dna2prot.sh \nwait\n')
			sf.write('source vfam.sh \nwait\n')
	if donrfilter:
		# if dorapsearch: sf.write(dirscr+'schedule2.py rap_nr.txt server.txt\nwait\n')
		# else: sf.write(dirscr+'schedule2.py blast_nr.txt server.txt\nwait\n')
		sf.write('source mergeSig.sh \nwait\n') #for diamond
		sf.write(dirscr+'schedule2.py diamond_nr.txt server.txt\nwait\n')
		#sf.write('source diamond_nr_view.sh \nwait\n')
		
	if doHmmer: 
		sf.write('source hmmer_nr_filter.sh >hmmernr.log  \nwait\n')
		sf.write('source hmmer_output_sort.sh >blastsort.log \nwait\n')
		
	else: 
		#sf.write('source blast_nr_filter.sh >blastnr.log  \nwait\n')
		sf.write('source diamond_nr_filter.sh >diamondnr.log  \nwait\n')
		sf.write('source blast_output_merge.sh \nwait\n')
		sf.write('source blast_output_sort.sh >blastsort.log \nwait\n')
	if doMyth: sf.write('source vfam_annot.sh\n')
	if genBlast: sf.write('source prepBlastFile.sh \n')
	sf.write('cat *.log > stats.logg\n')

	if doreAssemb: 
		sf.write('source reAssembleAll.sh \nwait\n')
		sf.write(dirscr+'firstpage.py reAssemble\nwait\n')
	else: sf.write(dirscr+'firstpage.py no\nwait\n')
	# if pair: sf.write('source prepPriceFile.sh \n')
	sf.write('echo '+password+' |find . "*" -print0 | sudo xargs -0 chmod 777\n')
	sf.write('source movetowww.sh \n')
	sf.write('source plot_pie.sh \n')
	sf.write('source plot_poly.sh  \nwait\n')
	sf.write('source mergeTable.sh \n')
	sf.write('echo '+password+' |find . "*" -print0 | sudo xargs -0 chmod 777\n')
	# sf.write('rsync -rP '+wd+'/'+base+' xdeng@dnasrv01dmzdr:/bsri/aspera_transfers/hold_test/\n')
	#if internal==False: sf.write('rsync -rP '+wd+'/'+base+' xdeng@dnasrv01dmzdr:/opt/aspera/shares/www/public/viral/\n')
	sf.close()
	cmd='echo '+password+' |find . "*" -print0 | sudo xargs -0 chmod 777'
	print cmd
	os.system(cmd)
