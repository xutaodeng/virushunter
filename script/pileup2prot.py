#!/usr/bin/env python
import sys

codontable = {
	'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
	'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
	'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
	'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
	'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
	'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
	'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
	'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
	'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
	'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
	'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
	'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
	'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
	'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
	'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
	'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
	}

def readFa(fafile):
	f=open(fafile, 'r')
	seq=[]
	for line in f:
		if not line.strip().startswith('>'):
			seq.append(line.strip())
	dna=''.join(seq)
	f.close()
	return dna.upper()

def pileup_dna2prot(inputfa, dnapile, protpile):
	dna=readFa(inputfa)
	f=open(dnapile,'r')
	of=open(protpile, 'w')
	header=f.readline().strip()
	of.write('\t'.join(header.split('\t')[0:3])+'\tprot\t'+'\t'.join(header.split('\t')[3:])+'\n')
	for line in f:
		parts = line.strip('\n').split('\t')
		id, pos, ref = parts[0:3]
		index = int(pos) #1 indexed
		frame = index%3
		calls=parts[3:]
		of.write(id+'\t'+pos+'\t'+ref+'\t')
		codonProt=''
		outs=[]
		for call in calls:
			out=[]
			syn = True #unsynomymous mutation
			for letter in call:
				if letter.upper() in set(['A','C','G','T']):
					if frame==0:
						codon=dna[index-3:index]
						mut=codon[0:2]+letter
					elif frame==1:
						codon=dna[index-1:index+2]
						mut=letter+codon[1:3]
					elif frame==2:
						codon=dna[index-2:index+1]
						mut=codon[0]+letter+codon[2]
					codonProt=codontable[codon]
					mutProt = codontable[mut]
					if mutProt != codonProt: syn=False
					out.append(mutProt)
				else: out.append(letter)
			if syn: out=[] #synomymous change
			outs.append(''.join(out))
		of.write(codonProt+'\t'+'\t'.join(outs)+'\n')
	of.close()
if __name__ == "__main__":
	pileup_dna2prot(sys.argv[1], sys.argv[2], sys.argv[3])