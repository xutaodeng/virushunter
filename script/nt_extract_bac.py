#!/usr/bin/env python

# wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz

#wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
#wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxcat.zip
# 7za e taxcat.zip
###############################################################################
# extract virus from nr and add to virus genome
###############################################################################
import gzip
from collections import defaultdict
import sys
import os

cat={}
f=open('categories.dmp', 'r')
for line in f:
	category, taxidspecies, taxid= line.strip().split()
	cat[taxid]=category
f.close()
print 'size of cat', len(cat)
f=gzip.open('nucl_gb.accession2taxid.gz', 'rb')
acc=set([])
i=0
f.readline()
for line in f:
	accession, accv, taxid, gi=line.strip().split()
	i+=1
	if cat.has_key(taxid) and cat[taxid]=='B':
		acc.add(accv)
	if i%100000000==0: print i, len(acc)
print 'size of accession', len(acc)
f.close()

# >X67214.1 B.taurus gene for adrenergic receptor beta 3
# >X13736.1 Bovine mRNA for adrenodoxin reductase C-term. (EC 1.18.1.2)
# >X05087.1 Cow Alu-type art2 sequence 7kb 5' of fetal globin gene
# >X05090.1 Cow alu-like art2 repetitive sequence 1kb 3' to fetal globin gene
# >X63475.1 B.taurus Alu-like repeat DNA
# >X63474.1 B.taurus Alu-like repeat DNA

def addTaxon():
	global acc
	f=gzip.open('nt.gz', 'rb')
	of=open('Bacteria.fa', 'w')

	bac=False
	b=0
	total=0
	for line in f:
		if line.strip().startswith('>'):
			total+=1
			if total%100000==0: print 'total fa sequence', total
			bac=False

			accession=line.strip().split()[0][1:]
			if accession in acc:
				bac=True
				print >>of, line.strip()
				b+=1
				if b%100000==0: print 'bac=', b
		elif bac: #sequence
			print >>of, line.strip()
	f.close()
	of.close()
	
##################################################################################################
# split Bacteria
##################################################################################################
def split_Bac():
    chunk=3200000000
    nchar=0
    chunkN=1
    filename = 'Bacteria'+str(chunkN)+'.fa'
    of = open(filename, 'w')
    print filename
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
                filename = 'Bacteria'+str(chunkN)+'.fa'
                of = open(filename, 'w')
                print filename
                print id, nchar
                print >>of, id.strip()
                print >> of, line.strip()
    of.close()
    f.close()

wd='/mnt/san/cluster/xdeng/nt/'
os.chdir(wd)
print 'current directory', os.getcwd()

print 'loading done'
addTaxon()
split_Bac()
#sys.exit()

ssh bsidna3 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria1.fa /mnt/cluster/xdeng/nt/Bacteria1 &
ssh bsidna4 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria2.fa /mnt/cluster/xdeng/nt/Bacteria2 &
ssh bsidna5 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria3.fa /mnt/cluster/xdeng/nt/Bacteria3 &
ssh bsidna6 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria4.fa /mnt/cluster/xdeng/nt/Bacteria4 &
ssh bsidna7 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria5.fa /mnt/cluster/xdeng/nt/Bacteria5 &
ssh bsidna8 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria6.fa /mnt/cluster/xdeng/nt/Bacteria6 &
ssh bsidna10 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria7.fa /mnt/cluster/xdeng/nt/Bacteria7 &
ssh bsidna11 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria8.fa /mnt/cluster/xdeng/nt/Bacteria8 &
ssh bsidna14 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria9.fa /mnt/cluster/xdeng/nt/Bacteria9 &
ssh bsidna15 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria10.fa /mnt/cluster/xdeng/nt/Bacteria10 &
wait
ssh bsidna3 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria11.fa /mnt/cluster/xdeng/nt/Bacteria11 &
ssh bsidna4 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria12.fa /mnt/cluster/xdeng/nt/Bacteria12 &
ssh bsidna5 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria13.fa /mnt/cluster/xdeng/nt/Bacteria13 &
ssh bsidna6 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria14.fa /mnt/cluster/xdeng/nt/Bacteria14 &
ssh bsidna7 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria15.fa /mnt/cluster/xdeng/nt/Bacteria15 &
ssh bsidna8 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria16.fa /mnt/cluster/xdeng/nt/Bacteria16 &
ssh bsidna10 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria17.fa /mnt/cluster/xdeng/nt/Bacteria17 &
ssh bsidna11 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria18.fa /mnt/cluster/xdeng/nt/Bacteria18 &
ssh bsidna14 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria19.fa /mnt/cluster/xdeng/nt/Bacteria19 &
ssh bsidna15 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria20.fa /mnt/cluster/xdeng/nt/Bacteria20 &
ssh bsidna19 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria21.fa /mnt/cluster/xdeng/nt/Bacteria21 &
ssh bsidna20 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria22.fa /mnt/cluster/xdeng/nt/Bacteria22 &
ssh bsidna23 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria23.fa /mnt/cluster/xdeng/nt/Bacteria23 &
ssh bsidna24 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria24.fa /mnt/cluster/xdeng/nt/Bacteria24 &
ssh bsidna25 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria25.fa /mnt/cluster/xdeng/nt/Bacteria25 &
ssh bsidna26 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria26.fa /mnt/cluster/xdeng/nt/Bacteria26 &
ssh bsidna30 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/nt/Bacteria27.fa /mnt/cluster/xdeng/nt/Bacteria27 &


#wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
#wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/mrna.fa.gz
# zcat hg38.fa.gz mrna.fa.gz > mrnadna.fa
# ssh bsidna12 /mnt/san/cluster/xdeng/tools/bowtie2-2.2.4/bowtie2-build /mnt/cluster/xdeng/hg38/mrnadna.fa /mnt/cluster/xdeng/hg38/mrnadna_bowtie &

