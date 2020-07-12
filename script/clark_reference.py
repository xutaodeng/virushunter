import gzip
from collections import defaultdict

f2=open('acc_tax_label.txt', 'r')
i=0
ca={}
for line in f2:
	i+=1
	if i%10000000==0: print i
	parts=line.strip().split()
	acc, cat, clas, fam, species =parts[1], parts[4], parts[5], parts[6], parts[7]
	ca[acc]='_'.join([cat, clas, fam, species])
f2.close()
print 'len acc', len(ca)

s=[]
# f1=open('samples.txt', 'r')
f1=open('nt.txt', 'r')
fs={}
counts=defaultdict(int)
ss=0
mapping=defaultdict(list)
for line in f1:
	fname=line.strip()
	ss+=1
	print ss, fname
	missing=0
	found=0
	fname='/mnt/cluster/xdeng/taxon/'+fname
	ff=gzip.open(fname, 'rb')
	good=False
	handle='_'
	for line in ff:
		if line.strip().startswith('>'):
			acc=line.strip().split()[0][1:]
			phage=False
			if 'PHAGE' in line.upper():
				phage=True
			try:
				cat=ca[acc]
				if phage:
					pa=cat.split('_')
					pa[0]='Phage'
					cat='_'.join(pa)
				fileNo=counts[cat]/100000
				handle=cat+'_'+str(fileNo)
				if counts[cat] %100000==0:
					print 'found', found, 'missing', missing
					fs[handle]=open('fasta/'+handle+'.fa', 'w')
					mapping[cat].append('/mnt/cluster/xdeng/taxon/fasta/'+handle+'.fa')
				fs[handle].write(line)
				counts[cat]+=1
				found+=1
				good=True
			except: 
				good=False
				missing+=1
		elif good:
			fs[handle].write(line)
	print 'found', found, 'missing', missing
	ff.close()
f1.close()

for key in fs.keys():
	fs[key].close()

of1=open('target1.txt', 'w')
of2=open('target2.txt', 'w')
of3=open('target3.txt', 'w')

handles=[of1, of2, of3]
hi=0
for key in mapping.keys():
	val=mapping[key]
	for v in val:
		of=handles[hi%3]
		hi+=1
		of.write(v+' '+key+'\n')
of1.close()
of2.close()
of3.close()
