#!/usr/bin/env python
import sys
import os
from collections import defaultdict
from operator import itemgetter, attrgetter
file1 = sys.argv[1]
label=sys.argv[2]
db = sys.argv[3]

f=open('test.tab', 'w')
f.close()
cmd='blastn -task megablast -num_threads 8 -max_target_seqs 1 -evalue 0.001 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen" '  +\
'-out test.tab -query '+file1+ ' -db /mnt/cluster/xdeng/blastdb/'+db #'vgenome'
os.system(cmd)
f=open('test.tab', 'r')
#'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'

def heights(hits):
	level=[]
	level_ends=len(hits)*[0]
	for i in xrange(len(hits)):
		hit = hits[i]
		start, end =hit[0], hit[1]
		for j in xrange(len(level_ends)):
			if start > level_ends[j]:
				level_ends[j] = end
				level.append(j)
				break
	return level

# def (hits):
	# for (s,e,l, qlen, alen) in hits:

def print_SVG():
	# header ='''
	# <!DOCTYPE html>
	# <html>
	# <body>
	# '''
	scale = 30
	for sid in vsvg.keys():
		maxh=maxhit(vsvg[sid])
		cov=coverage(vsvg[sid], 1000) 
		hits = sorted(vsvg[sid], key=itemgetter(0))
		slen = hits[0][2]
		ratio = slen/2000.0
		level = heights(hits)
		max_level=max(level)
		body=['<svg height="'+str(max_level*5+100)+'" width=2000">']
		body.append('<line x1="0" y1="0" x2="2000" y2="0" style="stroke:rgb(0,0,255);stroke-width:4" />')
		#body.append('<text x="0" y="5" fill="red">'+label+'_'+sid+'</text>')
		i=0
		len_align=0.0
		len_contig=0.0
		for (s,e,l, qlen, alen) in hits:
			if qlen >300 and alen > 200:
				len_contig+=qlen
				len_align+=min(alen, qlen) #alignment has indel and is greater than query
			y = str(10+level[i]*5)
			x1, x2=str(s/ratio), str(e/ratio)
			body.append('<line x1="'+x1+'" y1="'+y+'" x2="'+x2+'" y2="'+y+'" style="stroke:rgb(255,0,0);stroke-width:4" />')
			i+=1
		body.append('</svg><br>')
		chimera_Index='0'
		if len_contig>0: chimera_Index= str((len_contig-len_align)/len_contig)
		print label+' '+sid+' max_contig:'+str(maxh)+':chimera_index:'+chimera_Index+':c1000:'+str(cov)+'<br>'
	# trailer ='''</svg>
	# </body>
	# </html>
	# '''
		print '\n'.join(body)

#input a list of intervals, output coverage
def coverage(hits, thres):
	if len(hits)==0: return 0
	covset = set([])#cov = list(merge(hits))
	for (start, end,slen, qlen, alen) in hits:
		if alen < thres: continue
		for i in xrange(start, end):
			covset.add(i)
	return len(covset)

#input a list of intervals, output coverage
def maxhit(hits):
	if len(hits)==0: return 0
	covset = set([])#cov = list(merge(hits))
	max=0
	for (start, end,l, qlen, alen) in hits:
		if alen > max: max=alen
	return max

#vhits=defaultdict(list)
vsvg=defaultdict(list)
for line in f:
	parts = line.strip().split()
	qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, slen, qlen= parts[0:14] #length is alignment length
	sstart, send =  int(sstart), int(send)
	if sstart > send: sstart, send = send, sstart
	#vhits[sseqid].append((int(sstart), int(send)))
	#vhits[sseqid].append(int(length))#
	# chimera=False
	# if int(qlen) - int(length) >  100: 
		# chimera=True
	vsvg[sseqid].append((int(sstart), int(send), int(slen), int(qlen), int(length)))
f.close()

print_SVG()