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

def dna2prot(fafile, protfile, longest):
	f=open(fafile, 'r')
	of = open(protfile, 'w')
	seq=[]
	for line in f:
		if line.strip().startswith('>'):
			dna=''.join(seq)
			if dna!='':
				prots = translate_6(dna, longest)
				#print len(prots)
				i=1
				#print len(prots)
				for prot in prots:
					#print id, i, prot
					if len(prot) >=20:
						of.write(id+'prot_'+str(i)+'\n'+prot+'\n')
					i+=1
			seq=[]
			id=line.strip()
		else:
			seq.append(line.strip())
	dna=''.join(seq)
	prots = translate_6(dna, longest)
	#print len(prots)
	i=1
	for prot in prots:
		if len(prot) >=20:
			#print id, i, prot
			of.write(id+'prot_'+str(i)+'\n'+prot+'\n')
		i+=1
	f.close()
	of.close()

def revcomp(dna):
	rval=[]
	n=len(dna)
	for i in xrange(len(dna)):
		nuc = dna[n-i-1]
		if nuc=='A': rval.append('T')
		elif nuc=='C': rval.append('G')
		elif nuc=='G': rval.append('C')
		elif nuc=='T': rval.append('A')
	return ''.join(rval)

def translate(dna):
	n=len(dna)/3
	prot=[]
	for i in xrange(n):
		codon=dna[(i*3):(i*3+3)]
		try: aa=codontable[codon]
		except: return ''
		if aa=='_': break
		prot.append(aa)
	return ''.join(prot)

def translate_3(dna):
	return (translate(dna), translate(dna[1:]), translate(dna[2:]))
	
def translate_6(dna, longest):
	prot1, prot2, prot3 = translate_3(dna)
	prot4, prot5, prot6 = translate_3(revcomp(dna))
	prots = (prot1, prot2, prot3, prot4, prot5, prot6)
	if longest=="True":
		maxlen=-1
		maxprot=None
		for prot in prots:
			if len(prot) > maxlen:
				maxlen=len(prot)
				maxprot=prot
		return [maxprot]
	else: return prots
	
if __name__ == "__main__":
	try: longest=sys.argv[3]
	except: longest="False"
	dna2prot(sys.argv[1], sys.argv[2], longest)