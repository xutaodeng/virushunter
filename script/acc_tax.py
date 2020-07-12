#!/usr/bin/env python


###############################################################################
# extract virus from nr and add to virus genome
###############################################################################
import gzip
from collections import defaultdict
import sys
import os


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

	print 'annotations loaded'

def addTaxon():
	keep=set(['Archaea', 'Viruses', 'Viridiplantae', 'Fungi', 'Metazoa', 'Bacteria'])
	f=gzip.open('nucl_gb.accession2taxid.gz', 'rb')
	of= open('acc_tax_label.txt', 'w')
	total, good=0,0
	header=f.readline()
	print header
	for line in f:
		total+=1
		if total %1000000==0: print total, 'good', good
		parts=line.strip().split()
		acc, tid= parts[1], parts[2]
		try: path = nodes[tid].get_fullpath()
		except: path=[]
		category, family, clas, species= 'None', 'None','None','None'
		k=0
		for tid in path:
			k+=1
			try: level=nodes[tid].level.replace(' ', '_')
			except: level='None'
			try: nm=names[tid].replace(' ', '_')
			except: nm='None'
			#print tid, nm
			if level =='family':
				family=nm
			if level =='class':
				clas=nm
			if level =='species':
				species=nm
			if nm in keep:
				category=nm
				good+=1
		if category in keep:
			of.write(line.strip()+'\t'+category+'\t'+clas+'\t'+family+'\t'+species+'\n')

	print total, 'good=', good

loadTax()
addTaxon()