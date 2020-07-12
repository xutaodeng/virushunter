#!/usr/bin/env python

#0.	Add 7zip to path:
#a.	/mnt/cluster/tools/p7zip_9.20.1/bin

#1.	Back up the old nr and virus database to /mnt/cluster/xdeng/blastdb_date:
# cd /mnt/cluster/xdeng/blastdb run download.sh
# wget ftp://ftp.ncbi.nih.gov/refseq/release/viral/viral*protein.faa.gz
# wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz

# wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.zip
# wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_prot.zip
# wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxcat.zip
# wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip
# zcat viral*protein.faa.gz| gzip > viral.protein.fa.gz
# 7za e taxcat.zip
# 7za e gi_taxid_nucl.zip
# 7za e gi_taxid_prot.zip
# 7za e taxdmp.zip

#7.	Change mode
#a.	sudo chmod 777 * -R

###############################################################################
# extract virus from nr and add to virus genome
###############################################################################
import gzip
from collections import defaultdict
import sys
import os

wd='/mnt/san/cluster/xdeng/blastdb/'

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

giset=set([])
c1,c2, ap, h, v, nv, p, hm=0,0, 0, 0, 0, 0, 0, 0
cats=set([])
def addTaxon(infile, isnr):
	global giset, c1, c2, ap, h, v, nv, p, hm
	print 'isnr', isnr, 'giset size', len(giset)
	filter2={'Human immunodeficiency virus 1'.upper():0, 'Hepatitis C virus'.upper():0, 'Influenza A virus'.upper():0,'Hepatitis B virus'.upper():0}
	filter=filter2.keys()
	f=gzip.open(infile, 'rb')

	if isnr: 
		hf=open('human.virome.fa', 'a')
		virusf=open('virus.tmp.fa', 'a')
		phagef=open('phage.fa', 'a')
		nvf=open('nvnr.fa', 'a')
		df=open('diamond.fa', 'a')
	else:
		hf=open('human.virome.fa', 'w')
		virusf=open('virus.tmp.fa', 'w')
		phagef=open('phage.fa', 'w')
		nvf=open('nvnr.fa', 'w')
		df=open('diamond.fa', 'w')

	phage, herv, virus, appeared=False, False, False, False
	total=0
	for line in f:
		if line.strip().startswith('>'):
			total+=1
			if total%100000==0: print 'total fa sequence', total
			human, phage, herv, virus, appeared=False, False, False, False, False

			gi = line.strip()[4:].split('|')[0]
			if 'HUMAN' in line.upper() or 'SAPIENS' in line.upper():
				human=True
			if not isnr: #viral_protein
				giset.add(gi)
			else:  #nr
				for key in filter:
					if key in line.upper():
						filter2[key]+=1
						if filter2[key]> 1:
							ap+=1; appeared=True; break
				gs=line.strip().split('|')
				for g in gs:
					if g in giset: 
						ap+=1; appeared=True; break #This gi alread appeared

			if appeared: continue

			try: taxid = gis[gi]; c1+=1
			except: taxid ='0'
			try: path = nodes[taxid].get_fullpath()
			except: path=[]
			label=[]

			k=0
			category, family, genus, species= 'None', 'None','None','None'
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
					virus=True
				if 'PHAGE' in nm.upper():
					phage=True
				if 'HUMAN_ENDOGENOUS_RETROVIRUSES' == nm.upper():
					herv=True
			if herv:
				category='HERV'
			if phage:
				category='Phage'
			label.append('species'+'$'+species)
			label.append('genus'+'$'+genus)
			label.append('family'+'$'+family)
			label.append('category'+'$'+category)
			#label.append('subspecies'+'$'+subspecies)
			cats.add(category)
			if human and virus: hm+=1; print >>hf, line.strip()+' '+':'.join(label)
			if herv: h+=1; diamondLabel='HERV'; continue
			elif phage: diamondLabel='PHAGE'; p+=1; print >>phagef, line.strip()+' '+':'.join(label)
			elif virus: diamondLabel='VIRUS'; v+=1; print >>virusf, line.strip()+' '+':'.join(label)
			else: nv+=1; diamondLabel='NV'; print >>nvf, '>'+str(nv) #line.strip()+' '+':'.join(label)
			#print >>df, line.strip()+' '+':'.join(label)
			print >>df, '>'+diamondLabel+'_'+str(total)+'_'+':'.join(label) #line.strip()+' '+':'.join(label)
		else:#sequence
			if appeared: continue
			if human and virus:  print >>hf, line.strip().replace('-','')
			if herv: continue
			elif phage: print >>phagef, line.strip()
			elif virus: print >>virusf, line.strip().replace('-','')
			else: print >>nvf, line.strip()
			print >>df, line.strip().replace('-','')

	print cats
	c2=total-c1
	print 'duplicated', ap, 'with taxon', c1, 'without taxon', c2, 'virus', v, 'non-virus', nv, 'herv', h, 'phage', p, 'human', hm
	print filter2

def addLinlinHerv():
	f=open('virus.tmp.fa', 'r')
	f2=open('HERVaa.fasta', 'r')
	of=open('virus.fa', 'w')
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
print 'loading taxons'
loadTax()
print 'adding taxons virus'
addTaxon('viral.protein.fa.gz', False)
print 'adding taxons nr'
addTaxon('nr.gz', True)
addLinlinHerv()
os.system('segmasker -in virus.fa -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out virus_mask.asnb')
os.system('makeblastdb -in virus.fa -dbtype prot -parse_seqids -mask_data virus_mask.asnb -out virus_mask')
os.system('makeblastdb -in nvnr.fa -dbtype prot -parse_seqids -out nvnr')
os.system('segmasker -in phage.fa -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out phage_mask.asnb')
os.system('makeblastdb -in phage.fa -dbtype prot -parse_seqids -mask_data phage_mask.asnb -out phage_mask')
os.system('cat virus.fa phage.fa > virusphage.fa')
os.system('segmasker -in virusphage.fa -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out virusphage_mask.asnb')
os.system('makeblastdb -in virusphage.fa -dbtype prot -parse_seqids -mask_data phage_mask.asnb -out virusphage_mask')
os.system('diamond makedb --in diamond.fa -d diamond')
#add m13
#cat ./virus/virus.fa M13aa.fa > virus.tmp3.fa
#mv virus.tmp3.fa > ./virus/virus.fa