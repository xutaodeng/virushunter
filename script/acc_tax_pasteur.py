#!/usr/bin/env python


###############################################################################
# extract virus from nr and add to virus genome
###############################################################################
import gzip
from collections import defaultdict
import sys
import os
import re

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
			if child.level in ['subspecies', 'species', 'genus', 'family', 'order', 'superphylum', 'division', 'class', 'subphylum', 'phylum']: continue
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
accs={}
name2tid={}

def loadTax():
	global nodes, names, accs, name2tid
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
			name2tid[name]=tid
			names[tid]=name.replace(':', '_').replace('$', '_')
	f.close()
	of=open('tax_tree.txt', 'w')
	nodes['1'].printTree(0, of)
	of.close()

	print 'annotations loaded'

def addTaxon():
	f=open('rVDBv10.2_prot.fasta', 'r')
	virusf=open('viral_pasteur.fa', 'w')
	v=0
	goodtid=0
	nprot=0
	empty=False
	for line in f:
		if line.strip().startswith('>'):
			ss=[m.start() for m in re.finditer('[', line)]
			ee=[m.start() for m in re.finditer(']', line)]
			snames=
			#acc = line.strip().split('|')[0][1:]
			if line.strip().startswith('>|'): empty=True
			else: empty=False
			nprot+=1
			if nprot%100000==0: print 'nprot=', nprot, 'goodtid', goodtid
			parts=line.rsplit('[', 1)
			taxname=parts[1].strip()[0:-1]
			try: 
				tid=name2tid[taxname]
				goodtid+=1
			except: tid=None
			
			try: path = nodes[tid].get_fullpath()
			except: path=[]
		
			label=[]

			k=0
			category, family, genus, species= 'None', 'None','None','None'
			virus, phage=False, False
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

			if phage:
				category='Phage'
			label.append('species'+'$'+species)
			label.append('genus'+'$'+genus)
			label.append('family'+'$'+family)
			label.append('category'+'$'+category)
			#label.append('subspecies'+'$'+subspecies)
			#cats.add(category)
			# if herv: h+=1; diamondLabel='HERV'; continue
			#if phage: diamondLabel='PHAGE'; p+=1; print >>phagef, line.strip()+' '+':'.join(label)
			if virus and not phage and not empty: v+=1; print >>virusf, '>lcl|'+str(v)+'| '+line.strip()[1:]+' '+':'.join(label)
			# else: nv+=1; diamondLabel='NV'; print >>nvf, '>'+str(nv) #line.strip()+' '+':'.join(label)
			# #print >>df, line.strip()+' '+':'.join(label)
			# print >>df, '>'+diamondLabel+'_'+str(total)+'_'+':'.join(label) #line.strip()+' '+':'.join(label)
		else:#sequence
			if virus and not phage and not empty: print >>virusf, line.strip().replace('-','')

	print 'num virus', v
loadTax()
addTaxon()