###############################################################################
# extract virus from nr and add to virus genome
###############################################################################
import gzip
from collections import defaultdict
import sys

#################################################################################################
# step 1. run bowtie2 on the desired fastq or fasta file, e.g. caribu.fna
#################################################################################################
# ssh bsidna5 /mnt/cluster/tools/bowtie2-2.0.0-beta7/bowtie2 --no-hd --quiet --local -p 8 -x /mnt/cluster/xdeng/nt/nttax1 -f /mnt/san/cluster/xdeng/nt/caribu.fna -S /mnt/san/cluster/xdeng/nt/caribu_1.sam &
# ssh bsidna6 /mnt/cluster/tools/bowtie2-2.0.0-beta7/bowtie2 --no-hd --quiet --local -p 8 -x /mnt/cluster/xdeng/nt/nttax2 -f /mnt/san/cluster/xdeng/nt/caribu.fna -S /mnt/san/cluster/xdeng/nt/caribu_2.sam &
# ssh bsidna7 /mnt/cluster/tools/bowtie2-2.0.0-beta7/bowtie2 --no-hd --quiet --local -p 8 -x /mnt/cluster/xdeng/nt/nttax3 -f /mnt/san/cluster/xdeng/nt/caribu.fna -S /mnt/san/cluster/xdeng/nt/caribu_3.sam &
# ssh bsidna8 /mnt/cluster/tools/bowtie2-2.0.0-beta7/bowtie2 --no-hd --quiet --local -p 8 -x /mnt/cluster/xdeng/nt/nttax4 -f /mnt/san/cluster/xdeng/nt/caribu.fna -S /mnt/san/cluster/xdeng/nt/caribu_4.sam &
# ssh bsidna9 /mnt/cluster/tools/bowtie2-2.0.0-beta7/bowtie2 --no-hd --quiet --local -p 8 -x /mnt/cluster/xdeng/nt/nttax5 -f /mnt/san/cluster/xdeng/nt/caribu.fna -S /mnt/san/cluster/xdeng/nt/caribu_5.sam &
# ssh bsidna10 /mnt/cluster/tools/bowtie2-2.0.0-beta7/bowtie2 --no-hd --quiet --local -p 8 -x /mnt/cluster/xdeng/nt/nttax6 -f /mnt/san/cluster/xdeng/nt/caribu.fna -S /mnt/san/cluster/xdeng/nt/caribu_6.sam &
# ssh bsidna11 /mnt/cluster/tools/bowtie2-2.0.0-beta7/bowtie2 --no-hd --quiet --local -p 8 -x /mnt/cluster/xdeng/nt/nttax7 -f /mnt/san/cluster/xdeng/nt/caribu.fna -S /mnt/san/cluster/xdeng/nt/caribu_7.sam &
# ssh bsidna12 /mnt/cluster/tools/bowtie2-2.0.0-beta7/bowtie2 --no-hd --quiet --local -p 8 -x /mnt/cluster/xdeng/nt/nttax8 -f /mnt/san/cluster/xdeng/nt/caribu.fna -S /mnt/san/cluster/xdeng/nt/caribu_8.sam &
# ssh bsidna13 /mnt/cluster/tools/bowtie2-2.0.0-beta7/bowtie2 --no-hd --quiet --local -p 8 -x /mnt/cluster/xdeng/nt/nttax9 -f /mnt/san/cluster/xdeng/nt/caribu.fna -S /mnt/san/cluster/xdeng/nt/caribu_9.sam &
# ssh bsidna14 /mnt/cluster/tools/bowtie2-2.0.0-beta7/bowtie2 --no-hd --quiet --local -p 8 -x /mnt/cluster/xdeng/nt/nttax10 -f /mnt/san/cluster/xdeng/nt/caribu.fna -S /mnt/san/cluster/xdeng/nt/caribu_10.sam &
# ssh bsidna15 /mnt/cluster/tools/bowtie2-2.0.0-beta7/bowtie2 --no-hd --quiet --local -p 8 -x /mnt/cluster/xdeng/nt/nttax11 -f /mnt/san/cluster/xdeng/nt/caribu.fna -S /mnt/san/cluster/xdeng/nt/caribu_11.sam &
# ssh bsidna16 /mnt/cluster/tools/bowtie2-2.0.0-beta7/bowtie2 --no-hd --quiet --local -p 8 -x /mnt/cluster/xdeng/nt/nttax12 -f /mnt/san/cluster/xdeng/nt/caribu.fna -S /mnt/san/cluster/xdeng/nt/caribu_12.sam &


#################################################################################
#step 2. bowtie sam files only have gi, no tax
#load gi and species from nt_tax2 into dictionary tax
#Then run reads count: input sam files, output 
#################################################################################

f=open('nt_tax2.fa', 'r')
tax={}
print 'loading tax'
kk=0
for line in f:
	s=line.strip()
	if s.startswith('>'):
		s=s[1:]
		kk+=1
		if kk%1000000==0: print kk, s
		parts=s.split()
		tax[parts[0]]=s #gi is key, speci
f.close()

servers=['bsidna5','bsidna6','bsidna7','bsidna8','bsidna9','bsidna10',\
'bsidna11','bsidna12','bsidna13','bsidna14','bsidna15','bsidna16','bsidna17','bsidna18','bsidna19','bsidna20',\
'bsidna21','bsidna22','bsidna23','bsidna24','bsidna25','bsidna26','bsidna27','bsidna28','bsidna29', 'bsidna30']

nservers=len(servers)
sams =['caribu_1.sam','caribu_2.sam','caribu_3.sam','caribu_4.sam','caribu_5.sam','caribu_6.sam', \
'caribu_7.sam','caribu_8.sam','caribu_9.sam','caribu_10.sam','caribu_11.sam','caribu_12.sam' ]
counts =defaultdict(int)
ids=defaultdict()
for sam in sams:
	print sam
	f=open(sam, 'r')
	for line in f:
		parts = line.strip().split()
		readID, gi, cigar, read, score=parts[0], parts[2], parts[5], parts[9], parts[11]
		if gi=='*': continue
		species=tax[gi].split()[-1]
		score=int(score.split(':')[-1])
		if not ids.has_key(readID):
			ids[readID]=[species, score, tax[gi]]
		else:
			if score>ids[readID][-1]:
				ids[readID]=[species, score, tax[gi]]
	f.close()
for id in ids.keys():
	species=ids[id][0]
	counts[species]+=1
s=0
ofs={}
for (species, c) in counts.items():
	print species, c
	ofs[species]=open('caribou_'+species+'.fasta', 'w')
	s+=c
ofs['unaligned']=open('caribou_unaligned.fasta', 'w')
print 'total_aligned', s

f=open(sams[0], 'r')
ss=0
for line in f:
	parts = line.strip().split()
	ss+=1
	readID, species, cigar, read, score=parts[0], parts[2], parts[5], parts[9], parts[11]
	if ids.has_key(readID):
		species=ids[readID][0]
		labels=ids[readID][2]
		ofs[species].write( '>'+readID+'_'+labels+'\n'+read+'\n')
	else:
		ofs['unaligned'].write( '>'+readID+'_unaligned\n'+read+'\n')
f.close()
print 'total_reads', ss
for key in ofs.keys():
	ofs[key].close()

sys.exit()

##################################################################################################
# step 0-1: split annotated nt fastq file nt_tax2.fa into smaller chunks for bowtie2 index that can handle
# input nt_tax2.fa, output nttax1.fa, nttax2.fa ..., nttax12.fa
##################################################################################################
def split_NT():
    chunk=3200000000
    nchar=0
    chunkN=1
    #filename = 'nttax'+str(chunkN)+'.fa'
    filename = 'Bacteria'+str(chunkN)+'.fa'
    of = open(filename, 'w')
    print filename
    #f=gzip.open('nt_bowtie.fa.gz', 'rb')
    #f=open('nt_tax2.fa', 'r')
    f=open('Bacteria.fa', 'r')
    i=0
    for line in f:
        if line.strip().startswith('>'):
            i+=1
            if i%1000000==0: print i, line
            id=line
            print >>of, line.strip()
        else:
            print >>of, line.strip()
            nchar+=len(line.strip())
            if nchar > chunk:
                nchar=0
                of.close()
                chunkN+=1
                #filename = 'nttax'+str(chunkN)+'.fa'
                filename = 'Bacteria'+str(chunkN)+'.fa'
                of = open(filename, 'w')
                print filename
                print id, nchar
                print >>of, id.strip()
                print >> of, line.strip()
    of.close()
    f.close()

split_NT()
sys.exit()
############################################################################################
# step 0-2: make bowtie index on splitted nt
############################################################################################
# ssh bsidna3 bowtie2-build /mnt/cluster/xdeng/nt/Bacteria1.fa /mnt/cluster/xdeng/nt/Bacteria1 &
# ssh bsidna4 bowtie2-build /mnt/cluster/xdeng/nt/Bacteria2.fa /mnt/cluster/xdeng/nt/Bacteria2 &
# ssh bsidna5 bowtie2-build /mnt/cluster/xdeng/nt/Bacteria3.fa /mnt/cluster/xdeng/nt/Bacteria3 &
# ssh bsidna6 bowtie2-build /mnt/cluster/xdeng/nt/Bacteria4.fa /mnt/cluster/xdeng/nt/Bacteria4 &

# ssh bsidna3 bowtie2-build /mnt/cluster/xdeng/nt/nttax1.fa /mnt/cluster/xdeng/nt/nttax1 &
# ssh bsidna4 bowtie2-build /mnt/cluster/xdeng/nt/nttax2.fa /mnt/cluster/xdeng/nt/nttax2 &
# ssh bsidna5 bowtie2-build /mnt/cluster/xdeng/nt/nttax3.fa /mnt/cluster/xdeng/nt/nttax3 &
# ssh bsidna6 bowtie2-build /mnt/cluster/xdeng/nt/nttax4.fa /mnt/cluster/xdeng/nt/nttax4 &
# ssh bsidna7 bowtie2-build /mnt/cluster/xdeng/nt/nttax5.fa /mnt/cluster/xdeng/nt/nttax5 &
# ssh bsidna8 bowtie2-build /mnt/cluster/xdeng/nt/nttax6.fa /mnt/cluster/xdeng/nt/nttax6 &
# ssh bsidna9 bowtie2-build /mnt/cluster/xdeng/nt/nttax7.fa /mnt/cluster/xdeng/nt/nttax7 &
# ssh bsidna10 bowtie2-build /mnt/cluster/xdeng/nt/nttax8.fa /mnt/cluster/xdeng/nt/nttax8 &
# ssh bsidna11 bowtie2-build /mnt/cluster/xdeng/nt/nttax9.fa /mnt/cluster/xdeng/nt/nttax9 &
# ssh bsidna12 bowtie2-build /mnt/cluster/xdeng/nt/nttax10.fa /mnt/cluster/xdeng/nt/nttax10 &
# ssh bsidna13 bowtie2-build /mnt/cluster/xdeng/nt/nttax11.fa /mnt/cluster/xdeng/nt/nttax11 &
# ssh bsidna14 bowtie2-build /mnt/cluster/xdeng/nt/nttax12.fa /mnt/cluster/xdeng/nt/nttax12 &


###############################################################################
# step 0-0: adding taxons to the fasta: intput is nt.gz, output is nt_tax2.fa
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
            if nm in ['Bacteria']: continue
                
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
        names[tid]=name
f.close()
of=open('tax_tree.txt', 'w')
nodes['1'].printTree(0, of)
of.close()

gis={}
f=open('gi_taxid_nucl.dmp', 'r')
i=0
for line in f:
    i+=1
    if i%10000000==0: print 'loading gi_taxid', i
    gi, taxid = line.strip().split()
    if taxid=='0': continue
    gis[gi]=taxid
f.close()
print 'annotations loaded'

def filterNT(gis):
    f=gzip.open('nt.gz', 'rb')
    #of=gzip.open('nt_bowtie.fa.gz', 'wb')
    of=open('nt_tax2.fa', 'w')
    virus=set(['Viruses'])
    euk=set(['Viridiplantae','Fungi','Metazoa'])
    other_euk=set(['Eukaryota'])
    cell_viroid=set(['Bacteria','Archaea', 'Viroids'])
    fs={}
    i=0
    label=''
    for line in f:
        if line.strip().startswith('>'):
            i+=1
            if i%100000==0: print i
            gi = line.strip().split('|')[1]
            try: taxid = gis[gi]
            except: taxid ='0'
            try: path = nodes[taxid].get_fullpath()
            except: path=[]
            label='unclass'
            for tid in path:
                try: level=nodes[tid].level.replace(' ', '_')
                except: level='None'
                try: nm=names[tid].replace(' ', '_')
                except: nm='None'
                if nm in euk:
                    label=nm
                elif nm in other_euk:
                    label='other_euk'
                elif nm in cell_viroid:
                    label=nm
                elif nm in virus:
                    label=nm
                if label!='unclass':
                    break
            id_label=line.strip()[1:]+' '+label
            of.write('>'+id_label+'\n')
            if not fs.has_key(label):
               fs[label]=open(label+'.fa', 'w')
            fs[label].write('>'+id_label+'\n')
        else:
            of.write(line.strip()+'\n')
            fs[label].write(line.strip()+'\n')
    f.close()
    of.close()
    for label in fs.keys():
        fs[label].close()

filterNT(gis)
