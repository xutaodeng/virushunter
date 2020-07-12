###############################################################################
# extract virus from nr and add to virus genome
###############################################################################
import gzip
from collections import defaultdict
import sys

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
		

f=open('nodes.dmp', 'r')
nodes={}
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

names={}
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

gis={}
f=open('gi_taxid_nucl.dmp', 'r')
for line in f:
	gi, taxid = line.strip().split()
	if taxid=='0': continue
	gis[gi]=taxid
f.close()
print 'annotations loaded'

def addTaxonVirus(gis, f, of, of3):
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
				print >>of3, line.strip()+' '+':'.join(label)
		else:#sequence
			if found and not foundPhage and not foundHERV: 
				print >>of, line.strip()
			else:
				print >>of3, line.strip()
			# if foundHERV or foundPhage: 
				# print >>of2, line.strip()
	print cats
	print 'with name', c1, 'without name', c2, 'virus', v, 'non-virus', nv, 'herv', h, 'phage', p

f2=gzip.open('viral.1.1.genomic.fna.gz', 'rb')
of=open('virus_dna.fa', 'w')
of3=open('virus_dna_noLabel.fa', 'w')
print 'adding taxons virus'
addTaxonVirus(gis, f2, of, of3)
f2.close()
of.close()
of3.close()