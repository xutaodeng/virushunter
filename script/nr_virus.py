#!/usr/bin/env python

#0.	Add 7zip to path:
#a.	/mnt/cluster/tools/p7zip_9.20.1/bin



#1.	Back up the old nr and virus database:

#mv /mnt/san/cluster/xdeng/blastdb/nr /mnt/san/cluster/xdeng/blastdb/nr_today
#mv /mnt/san/cluster/xdeng/blastdb/virus /mnt/san/cluster/xdeng/blastdb/virus_today

#2.	Download the virus genome and nr to current directory /mnt/san/cluster/xdeng/blastdb/

# wget ftp://ftp.ncbi.nih.gov/refseq/release/viral/viral.1.protein.faa.gz
# wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
## wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/microbial/*genomic.fna.gz 

# cd refseq_download
# wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/*nonredundant*protein.faa.gz -P ./refseq_download
# zcat *nonredundant*protein.faa.gz > refseq.fa

# #3.	Update taxonomy:
# wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.zip
# wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_prot.zip
# wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxcat.zip
# wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip
# 7za e taxcat.zip
# 7za e gi_taxid_nucl.zip
# 7za e gi_taxid_prot.zip
# 7za e taxdmp.zip
## 7za e nr.zip
# #4.	Create new empty directory
# mkdir nr
# mkdir virus
# #5.	Process the downloaded files 
                # python nr_virus.py
# #6.	Run makeblastdb
# segmasker -in virus.fa -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out virus_mask.asnb
# makeblastdb -in virus.fa -dbtype prot -parse_seqids -mask_data virus_mask.asnb -out virus_mask
# makeblastdb -in nvrefseq.fa -dbtype prot -parse_seqids -out nvrefseq


# # segmasker -in phageHERV.fa -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out phageHERV_mask.asnb
# # makeblastdb -in phageHERV.fa -dbtype prot -parse_seqids -mask_data phageHERV_mask.asnb -out phageHERV_mask
# # segmasker -in nvrefseq.fa -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out nvrefseq_mask.asnb
# # makeblastdb -in nvrefseq.fa -dbtype prot -parse_seqids -mask_data nvrefseq_mask.asnb  -out nvrefseq_mask

 # prerapsearch -d virus.fa -n virusrap
 # prerapsearch -d nvrefseq.fa -n nvrap

#7.	Change mode
#a.	sudo chmod 777 * -R



###############################################################################
# extract virus from nr and add to virus genome
###############################################################################
import gzip
from collections import defaultdict
import sys
import os

dirscr='/mnt/cluster2/xdeng/script/'
wd='/mnt/san/cluster/xdeng/blastdb/'
def virusdb():
	nrfilter=['Human immunodeficiency virus 1'.upper(), 'Hepatitis C virus'.upper(), 'Influenza A virus'.upper(),'Hepatitis B virus'.upper()]
	vf=gzip.open('viral.1.protein.faa.gz', 'rb')
	nf=gzip.open('nr.gz', 'rb')
	id =''

	#of2=open('refseq/nvrefseq.fa','w')
	of=open('virus.tmp.fa', 'w')
	virus, skip=False, False
	c,x,y=0,0,0
	gis0=set()
	for line in vf:
		if line.strip().startswith('>'):
			skip=False
			gi=line.strip().split('|')[1]
			gis0.add(gi)
			# if 'PHAGE' in line.upper():
				# skip=True
			# else:
			c+=1
			print >>of, line.strip()
		else:
			if skip: continue
			else: print >>of, line.strip()
	print 'virus loaded', 'virusdb', c
	k=0
	z=0
	j=0
	for line in nf:
		if line.strip().startswith('>'):
			j+=1
			#if j%1000000==0: print j
			id=line.strip().upper()
			gilist=line.strip().split('|')
			skip=False
			for gi in gilist:
				if gi in gis0:
					k+=1
					#print k, gi
					skip=True
					break
			if ('VIRUS' in id or 'VIRAL' in id):
				virus=True
				for junk in nrfilter:
					if junk in id:
						skip=True
						z+=1
						break
			else: virus=False
			if not skip and virus:
				x+=1
				print >>of, line.strip()
			elif not skip:
				y+=1
				#print >>of2, line.strip()
		elif not skip and virus:
			print >>of, line.strip()
		elif not skip:
			pass
			#print >>of2, line.strip()
	vf.close()
	nf.close()
	of.close()
	#of2.close()
	print 'virusdb', c, 'nvirus', x, 'totalvirus', c+x, 'non-virus', y, 'skipped_overlap', k, 'nr_filter', z
	
###############################################################################
# adding taxons to the fasta
###############################################################################
class Node(object): 
	def __init__(self, tid, level=None): 
		self.tid = tid 
		self.parent = None
		self.children = []
		self.level=level

	def add_parent(self, parent):
		self.parent = parent
		parent.children.append(self)
		
	def get_fullpath(self):
		rval=[self.tid]
		node=self
		while node.parent!=None:
			tmp=node.parent
			rval.append(tmp.tid)
			node=tmp
		return rval
	
	def printTree(self, depth, of):
		for child in self.children:
			# if child.level in ['subspecies', 'species', 'genus', 'family', 'order', 'superphylum', 'division', 'class', 'subphylum', 'phylum']: continue
			level=child.level
			if level=='no rank': level=''
			try:
				nm=names[child.tid].replace(' ', '_')
			except:
				nm='none'
			print >>of, "\t" * depth, '-', depth, level, nm
			#if nm in ['Viruses']: continue
				
			child.printTree(depth+1, of)

nodes={}
names={}
gis={}
def loadTax():
	global nodes, names, gis
	f=open('nodes.dmp', 'r')
	
	for line in f: #create nodes
		parts = line.strip().split('|')
		tid, pid, level = parts[0].strip(), parts[1].strip(), parts[2].strip()
		nodes[tid]= Node(tid, level)
	f.close()

	f=open('nodes.dmp', 'r')
	for line in f:
		parts = line.strip().split('|')
		tid, pid, level = parts[0].strip(),parts[1].strip(),parts[2].strip()
		n1 = nodes[tid]
		if tid!=pid:
			n2 = nodes[pid]
			n1.add_parent(n2)
	f.close()

	f=open('names.dmp', 'r')
	for line in f:
		if 'scientific name' in line:
			parts=line.strip().split('|')
			tid=parts[0].strip()
			name=parts[1].strip()
			names[tid]=name.replace(':', '_').replace('$', '_')
	f.close()
	of=open('tax_tree.txt', 'w')
	nodes['1'].printTree(0, of)
	of.close()

	f=open('gi_taxid_prot.dmp', 'r')
	for line in f:
		gi, taxid = line.strip().split()
		if taxid=='0': continue
		gis[gi]=taxid
	f.close()
	print 'annotations loaded'

def addTaxonVirus(f, of, of2):
	c1,c2=0,0
	h, v, nv, p = 0, 0, 0, 0
	cats=set([])
	for line in f:
		if line.strip().startswith('>'):
			gi = line.strip().split('|')[1]
			try: taxid = gis[gi]; c1+=1
			except: taxid ='0'; c2+=1
			try: path = nodes[taxid].get_fullpath()
			except: path=[]
			label=[]

			k=0
			found=False
			foundHERV=False
			foundPhage=False
			category, family, genus, species= \
					  'None', 'None','None','None'
			for tid in path:
				k+=1
				try: level=nodes[tid].level.replace(' ', '_')
				except: level='None'
				try: nm=names[tid].replace(' ', '_')
				except: nm='None'
				
				if level =='species':
					species=nm
				elif level =='genus':
					genus=nm
				elif level =='family':
					family=nm
				elif level=='no_rank' and k ==(len(path)-2):
					category=nm
				elif level =='superkingdom' and nm=='Viruses':
					found=True
				if 'PHAGE' in nm.upper():
					foundPhage=True
				if 'Human_endogenous_retroviruses' == nm:
					foundHERV=True
			if foundHERV:
				category='HERV'
			if foundPhage:
				category='Phage'
			label.append('species'+'$'+species)
			label.append('genus'+'$'+genus)
			label.append('family'+'$'+family)
			label.append('category'+'$'+category)
			#label.append('subspecies'+'$'+subspecies)
			cats.add(category)
			# if foundHERV or foundPhage:
				# if foundHERV: h+=1
				# if foundPhage: p+=1
				# print >>of2, line.strip()+' '+':'.join(label)
			if found and not foundPhage and not foundHERV:
				v+=1 #virus
				print >>of, line.strip()+' '+':'.join(label)
			else:
				nv+=1 #non-vurs
				#print >>of3, line.strip()+' '+':'.join(label)
		else:#sequence
			if found and not foundPhage and not foundHERV: 
				print >>of, line.strip()
			#else:
			#	print >>of3, line.strip()
			# if foundHERV or foundPhage: 
				# print >>of2, line.strip()
	print cats
	print 'with name', c1, 'without name', c2, 'virus', v, 'non-virus', nv, 'herv', h, 'phage', p

def addTaxonNV(f, of3):
	c1,c2=0,0
	v, nv = 0, 0
	cats=set([])
	for line in f:
		if line.strip().startswith('>'):
			gi = line.strip().split('|')[1]
			try: taxid = gis[gi]; c1+=1; annot=True
			except: taxid ='0'; c2+=1; annot=False
			try: path = nodes[taxid].get_fullpath()
			except: path=[]
			label=[]

			k=0
			found=False
			category, family, genus, species= 'None', 'None','None','None'
			for tid in path:
				k+=1
				try: level=nodes[tid].level.replace(' ', '_')
				except: level='None'
				try: nm=names[tid].replace(' ', '_')
				except: nm='None'
				#label2.append(level+'$'+nm)
				
				#if level =='subspecies':
				 #   subspecies=nm
				if level =='species':
					species=nm
				elif level =='genus':
					genus=nm
				elif level =='family':
					family=nm
				elif level=='no_rank' and k ==(len(path)-2):
					category=nm
				elif level =='superkingdom' and nm=='Viruses':
					found=True

			label.append('species'+'$'+species)
			label.append('genus'+'$'+genus)
			label.append('family'+'$'+family)
			label.append('category'+'$'+category)
			if found:
				v+=1 #virus
				#print >>of, line.strip()+' '+':'.join(label)
			if annot and not found:
				nv+=1 #non-vurs
				print >>of3, '>'+str(nv)#line.strip()+' '+':'.join(label)
		else:#sequence
			if annot and not found:
				print >>of3, line.strip()
	print 'with name', c1, 'without name', c2, 'virus', v, 'non-virus', nv
def addLinlinHerv():
	f=open('virus.tmp2.fa', 'r')
	f2=open('HERVaa.fasta', 'r')
	try: os.mkdir('virus')
	except: pass
	of=open('virus/virus.fa', 'w')
	for line in f:
		of.write(line)
	for line in f2:
		if line.strip()=='':
			continue
		elif line.strip().startswith('>'):
			if line.strip()[1:].strip().startswith('gi|'):
				species =  '_'.join(line.strip()[1:].split()[1:])
			else:
				species='_'.join(line.strip()[1:].strip().split())
			of.write(line.strip()+' species$'+species+':genus$HERV:family$HERV:category$HERV\n')
				#species$Lettuce_necrotic_stunt_virus:genus$Tombusvirus:family$Tombusviridae:category$ssRNA_viruses
		else:
			of.write(line)
	f.close()
	f2.close()
	of.close()

os.chdir(wd)
print 'current directory', os.getcwd()
#virusdb()
#sys.exit()

print 'loading taxons'
loadTax()

print 'adding taxons virus'
f2=open('virus.tmp.fa', 'r')
of=open('virus.tmp2.fa', 'w')
of2=open('phageHERV.fa', 'w')

addTaxonVirus(f2, of, of2)
f2.close()
of2.close()
of.close()
f=open('refseq_download/refseq.fa', 'r')
of3=open('nvrefseq.fa', 'w')
print 'adding taxons NV'
addTaxonNV(f, of3)
of3.close()
addLinlinHerv()
#add m13
#cat ./virus/virus.fa M13aa.fa > virus.tmp3.fa
#mv virus.tmp3.fa > ./virus/virus.fa