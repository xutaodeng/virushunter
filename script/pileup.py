#!/usr/bin/env python

import sys
from collections import defaultdict

standardCodon={
"TTT": ("F", "Phe"),	"TTC": ("F", "Phe"),	"TTA": ("L", "Leu"),	"TTG": ("L", "Leu"),	
"TCT": ("S", "Ser"),	"TCC": ("S", "Ser"),	"TCA": ("S", "Ser"),	"TCG": ("S", "Ser"),	
"TAT": ("Y", "Tyr"),	"TAC": ("Y", "Tyr"),	"TAA": ("*", "Ter"),	"TAG": ("*", "Ter"),	
"TGT": ("C", "Cys"),	"TGC": ("C", "Cys"),	"TGA": ("*", "Ter"),	"TGG": ("W", "Trp"),	

"CTT": ("L", "Leu"),	"CTC": ("L", "Leu"),	"CTA": ("L", "Leu"),	"CTG": ("L", "Leu"),	
"CCT": ("P", "Pro"),	"CCC": ("P", "Pro"),	"CCA": ("P", "Pro"),	"CCG": ("P", "Pro"),	
"CAT": ("H", "His"),	"CAC": ("H", "His"),	"CAA": ("Q", "Gln"),	"CAG": ("Q", "Gln"),	
"CGT": ("R", "Arg"),	"CGC": ("R", "Arg"),	"CGA": ("R", "Arg"),	"CGG": ("R", "Arg"),	

"ATT": ("I", "Ile"),	"ATC": ("I", "Ile"),	"ATA": ("I", "Ile"),	"ATG": ("M", "Met"),	
"ACT": ("T", "Thr"),	"ACC": ("T", "Thr"),	"ACA": ("T", "Thr"),	"ACG": ("T", "Thr"),	
"AAT": ("N", "Asn"),	"AAC": ("N", "Asn"),	"AAA": ("K", "Lys"),	"AAG": ("K", "Lys"),	
"AGT": ("S", "Ser"),	"AGC": ("S", "Ser"),	"AGA": ("R", "Arg"),	"AGG": ("R", "Arg"),	

"GTT": ("V", "Val"),	"GTC": ("V", "Val"),	"GTA": ("V", "Val"),	"GTG": ("V", "Val"),	
"GCT": ("A", "Ala"),	"GCC": ("A", "Ala"),	"GCA": ("A", "Ala"),	"GCG": ("A", "Ala"),	
"GAT": ("D", "Asp"),	"GAC": ("D", "Asp"),	"GAA": ("E", "Glu"),	"GAG": ("E", "Glu"),	
"GGT": ("G", "Gly"),	"GGC": ("G", "Gly"),	"GGA": ("G", "Gly"),	"GGG": ("G", "Gly"),	
}

ff=open(sys.argv[5], 'r')
header=ff.readline()
refseq=['N'] #1 indexed
for line in ff:
	if line.strip().startswith('>'): continue
	refseq.append(line.strip())
refseq=''.join(refseq)

f=open(sys.argv[1], 'r') #pileup file
f1=open(sys.argv[2], 'r') #script file
of1=open(sys.argv[3],'w')
of2=open(sys.argv[4],'w')
of3=open('consensus.fasta', 'w')

samples=f1.readline().split()[7:-2]
samples=[s.split('.')[0] for s in samples]
print samples
n = len(samples)
print 'nsample', n
seqs=[]
for i in xrange(n): 
	#seqs.append(list(refseq))
	seqs.append(['N']*len(refseq))

k=0
of1.write("genome\tpos\tref\t"+'\t'.join(samples)+'\n')
of2.write("genome\tpos\tref\t"+'\t'.join(samples)+'\n')
for line in f:
	k+=1
	if k%100==0: print 'line', k
	# if k<90: continue
	# if k==100:  break; #sys.exit();continue
	parts=line.strip('\n').split('\t')
	chro, pos, ref = parts[0:3]
	pos=int(pos)
	data = parts[3:]
	of1.write('\t'.join(parts[0:3]))
	of2.write('\t'.join(parts[0:3]))
	for i in xrange(n):
		coverage, read, qual = data[i*3:(i*3+3)]
		coverage = int(coverage)
		read =read.upper()
		same, A, C, G, T, N, plus, minus,star= 0, 0,0,0,0,0,0,0,0
		plusnumbers, minusnumbers=defaultdict(int),defaultdict(int)
		j=0
		while j < len(read):
			r=read[j]
			if r=='^': j+=2
			elif r=='.' or r==',': same+=1; j+=1
			elif r=='+' or r=='-': 
				if not read[j+1].isdigit(): j+=1; continue
				if r=='+': plus+=1
				else: minus+=1
				j+=1
				number=[]
				while read[j].isdigit():
					number.append(read[j])
					j+=1
				try: number=int(''.join(number))
				except: print j, r, number, read; sys.exit(1)
				sub=read[j:(j+number)].upper()
				if r=='+': plusnumbers[sub]+=1; # print 'sub',sub, 'number', number, read; 
				else: minusnumbers[sub]+=1
				j+=number
			elif r=='A': A+=1; j+=1
			elif r=='C': C+=1; j+=1
			elif r=='G': G+=1; j+=1
			elif r=='T': T+=1; j+=1
			elif r=='N': N+=1; j+=1
			elif r=='*': star+=1; j+=1
			else: j+=1
		mm1, subplus=0,None
		for key, value in plusnumbers.items():
			if value > mm1:
				mm1=value
				subplus=key
		mm2, subminus=0,None
		for key, value in minusnumbers.items():
			if value > mm2:
				mm2=value
				subminus=key

		cons =max([same, A, C, G, T, N, star, mm1, mm2])
		if coverage==0: alt='N'
		elif same==cons:alt=ref
		elif A==cons: alt='A'
		elif C==cons: alt='C'
		elif G==cons: alt='G'
		elif T==cons: alt='T'
		elif N==cons: alt='N'#; print 'lineno', k, read, ref, same, A, C, G, T, N, "alt=",alt,  coverage, mm1, subplus, plusnumbers, minusnumbers; #sys.exit()
		elif mm1==cons: alt=ref+subplus 
		elif mm2==cons: alt=''
		elif star==cons: alt=''#; print read
		
		if coverage>0 and mm1 > 0.8*coverage:
			alt=ref+subplus 
		#if len(plusnumbers)>0: print ref, same, A, C, G, T,alt,  coverage, read[0:30], mm1, subplus, plusnumbers;
		#print 'line', k , read, ref, same, A, C, G, T, N, "alt=",alt,  coverage, mm1, subplus, plusnumbers, minusnumbers; sys.exit()
		if len(alt)>1:
			 print 'line', k , ref, 'same=',same, 'A=',A, 'C=',C, 'G=',G, 'T=',T, 'N=',N, '*=',star, "alt=",alt,  'coverage=',coverage, mm1, subplus, plusnumbers, minusnumbers, read[0:30]
		seqs[i][pos]=alt
		# if len(alt)>1 or alt=='' or alt=='N':  
			# if len(alt)>1: print 'Insertion'
			# elif alt=='N': 
				# print "N"
				# print 'line', k , ref, 'same=',same, 'A=',A, 'C=',C, 'G=',G, 'T=',T, 'N=',N, '*=',star, "alt=",alt,  coverage, mm1, subplus, plusnumbers, minusnumbers, read[0:30]
				# print i, pos, seqs[i][pos]
		rval, rval2=[],[]
		cthres, freqthres=0,0 #0.005
		keep=False
		if coverage >cthres:
			if A/coverage>freqthres: 
				rval.append(str(A)+'A')
				rval2.append(str(round(A/coverage*100,1))+'A')
				keep=True
			if C/coverage>freqthres: 
				rval.append(str(C)+'C')
				rval2.append(str(round(C/coverage*100,1))+'C')
				keep=True
			if G/coverage>freqthres: 
				rval.append(str(G)+'G')
				rval2.append(str(round(G/coverage*100,1))+'G')
				keep=True
			if T/coverage>freqthres: 
				rval.append(str(T)+'T')
				rval2.append(str(round(T/coverage*100,1))+'T')
				keep=True
			if plus/coverage>freqthres: 
				rval.append(str(plus)+'+')
				rval2.append(str(round(plus/coverage*100,1))+'+')
				keep=True
			if minus/coverage>freqthres: 
				rval.append(str(minus)+'-')
				rval2.append(str(round(minus/coverage*100,1))+'-')
				keep=True
			if same/coverage>freqthres:
				rval.append(str(same)+ref)
				rval2.append(str(round(same/coverage*100,1))+ref)
		rval=''.join(rval)
		rval2=''.join(rval2)
		if keep: of1.write('\t'+rval)
		else: of1.write('\t')
		if keep: of2.write('\t'+rval2)
		else: of2.write('\t')
	of1.write('\n')
	of2.write('\n')

for i in xrange(len(seqs)):
	of3.write('>'+samples[i]+'\n')
	j=1
	lenseq=0
	while j  < len(seqs[i]):
		of3.write(seqs[i][j])
		lenseq+=len(seqs[i][j])
		#if i==101 and j==101: print seqs[i][j]
		if len(seqs[i][j])>1: print 'insertion', seqs[i][j]
		#elif len(seqs[i][j])=='N': pass; print 'N', 'i=',i, 'j=',j
		j+=1
	print i, samples[i], 'lenseq', lenseq
	of3.write('\n')

of3.close()
f.close()
f1.close()
of1.close()
of2.close()
