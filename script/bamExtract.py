muts=[ (17,62006799), (6,29911916), (3,51978220), (3,158370034), (3,184037533), (19,50168927)]
window=100
# for line in f:
	# muts.append(line.strip().split())
bams=['Pt1', 'Pt2',  'Xeno2', 'Xeno31', 'Xeno32', 'Xeno4', 'Xeno5', 'Xeno6', 'Xeno7' ]

for (chro, pos) in muts:
	pos=int(pos)
	chro=str(chro)
	for bam in bams:
		cmd='samtools view -h ../'+bam+'/2pass/s6.bam '+chro+':'+str(pos-window)+'-'+str(pos+window)+' > '+bam+'_chr'+chro+'_'+str(pos)+'.sam'
		print cmd