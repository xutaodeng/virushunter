#!/usr/bin/env python
import sys
blast = {}
f=open(sys.argv[1], 'r')
for line in f:
	if 'max_contig' in line:
		label, ref, measure = line.strip().strip('<br>').split()
		if label.count('_')>4: #test Kmer
			prog, kmer, dat = label.split(':')[1].split('_')
			if '_' in dat: dat=dat.split('_')[1]
			max_con=measure.split(':')[1]
			c1000=measure.split(':')[5]
			chimera_index=measure.split(':')[3]
			blast[dat+' '+prog+' '+kmer+' '+ref] = (max_con, chimera_index, c1000)
		else: #not test kmer
			label1 = label.split('/')[0].split('_')[1]
			#print label
			prog, dat = label.split(':')[1].split('_',1)
			if '_' in dat: dat=dat.split('_')[1]
			if label1 !='individual':
				prog = label1
			if prog=='soap': prog = 'S'
			elif prog=='meta': prog='V'
			elif prog=='abyss': prog = 'A'
			elif prog=='mira': prog='M'
			elif prog=='omega': prog='G'
			elif prog=='celera': prog='W'
			elif prog=='trinity': prog='T'
			max_con=measure.split(':')[1]
			chimera_index=measure.split(':')[3]
			c1000=measure.split(':')[5]
			blast[dat+' '+prog+' '+ref] = (max_con, chimera_index, c1000)
f.close()

for id in blast.keys():
	print id, blast[id][0], blast[id][1], blast[id][2]

