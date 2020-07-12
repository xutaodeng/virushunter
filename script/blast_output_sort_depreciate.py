#!/usr/bin/env python

from collections import defaultdict
import operator
import sys
import os
import os.path
import linecache
import re, string
def outputHeader(of):
	head ='''
	<html>
<script src="../../sorttable.js">
</script>
<script src="../../ajax_select.js">
</script>
<input type="checkbox" name='checkall' value='checkall' onclick="checkAll(this)">
Check All <button onclick="exportData()" >Export Selected</button> <br>
 <div id="txtHint"> </div>
<table class="sortable">
<tr><th class="sorttable_nosort"></th><th>Category</th><th>Family</th><th>Genus</th><th>Species</th><th>Virus</th><th class="sorttable_numeric">Low_V_Evalue</th><th class="sorttable_numeric">Low_NR_Evalue</th><th class="sorttable_numeric">Hits</th><th class="sorttable_nosort"></th></tr>
'''
#<tr><th class="sorttable_nosort"><th>Virus</th><th class="sorttable_numeric">Low_V_Evalue</th><th class="sorttable_numeric">Low_NR_Evalue</th><th class="sorttable_numeric">Hits</th><th class="sorttable_nosort"></th>


	of.write(head+'\n')

def CacheLines(input): 
	cache=defaultdict(list)
	VE={}
	NE={}
	tax={}
	f = open(input, 'r')
	i=0
	for line in f:
		i+=1
		if i%11 == 4:
			taxonomy=line.strip().split()[-1]
			index=line.rfind('[')
			if index==-1:
				parts=line.upper().split()
				for k in xrange(len(parts)):
					part=parts[k]
					if part=='VIRUS' or part=='VIRUS,' or part=='VIRUS.':
						part=part.strip('.').strip(',')
						try: virus='['+parts[k-1]+' '+part+']'
						except: print parts; virus='['+part+']'
						break
					elif part.find('VIRUS') !=-1:
						virus='['+part+']'
						break
			else:
				index2=line.rfind(']')
				virus=line.strip()[index:(index2+1)].upper()
				parts = virus.split()
				try: 
					if parts[-2]=='VIRUS':
						virus=' '.join(parts[0:-1])+']'
				except:
					pass
			cache[virus].append(i-3)
			tax[virus]=taxonomy
		elif i%11==6: #veval
			ve=float(line.split()[-1])
			if VE.has_key(virus) and ve<VE[virus]:
				VE[virus]=ve
			elif not VE.has_key(virus):
				VE[virus]=ve
		elif i%11==7: #nr eval, LVNE
			ne=line.split()[-1]
			try: ne=float(ne)
			except: ne=1
			if NE.has_key(virus) and ne<NE[virus]:
				NE[virus]=ne
			elif not NE.has_key(virus):
				NE[virus]=ne
	f.close()
	return cache, VE, NE, tax

def printVirusSummary(cache, VE, NE, tax, input):
	of = open(os.path.splitext(input)[0]+'.html', 'w')
	keys=sorted(cache.keys())
	outputHeader(of)
	numvirus=0
	for virus in keys:
		try: 
			species, genus, family, cat = tax[virus].split(':')[0:4]
			cat=cat.split('$')[1]
			family=family.split('$')[1]
			genus=genus.split('$')[1]
			species=species.split('$')[1]
		except: cat,species, genus, taxonomy = 'NA', 'NA', 'NA', 'NA'
		base=os.path.basename(input)
		virname='_'.join(virus.strip('[').strip(']').replace('/','_').split())[0:40]
		virname = re.sub('[\W_]+', '_', virname)
		label = base+'_'+virname
		f1 = 'aln/'+label+'.txt'
		f2 = 'fasta/'+label+'.fa'
		of.write('<tr><td><input type="checkbox" value="'+label+'"> </td>')
		#of.write('<td>'+virname+'</td><td>'+str(VE[virus])+'</td><td>'+str(NE[virus])+'</td><td><a href=\"'+f1+'\">'+str(len(cache[virus]))+'</a></td><td>')
		of.write('<td>'+cat+'</td><td>'+family+'</td><td>'+genus+'</td><td>'+species+'</td><td>'+virname+'</td><td>'+str(VE[virus])+'</td><td>'+str(NE[virus])+'</td><td><a href=\"'+f1+'\">'+str(len(cache[virus]))+'</a></td><td>')
		of.write('<a href=\"'+f2+'\">fasta</a></td></tr>\n')
		numvirus+=1
	of.write('</table></html>\n')
	of.close()
	print input, 'num_virus =', numvirus

def OutputVirus(cache, input):
	keys=sorted(cache.keys())
	try: os.mkdir(os.path.dirname(input)+'/aln/')
	#print 'dir', os.path.dirname(input)+'/aln/'
	except: pass
	try: os.mkdir(os.path.dirname(input)+'/fasta/')
	except: pass
	for virus in keys:
		ids=set()
		base=os.path.basename(input)
		virname='_'.join(virus.strip('[').strip(']').replace('/','_').split())[0:40]
		virname = re.sub('[\W_]+', '_', virname)
		of = open(os.path.dirname(input)+'/aln/'+base+'_'+virname+'.txt', 'w')
		of2 = open(os.path.dirname(input)+'/fasta/'+base+'_'+virname+'.fa', 'w')
		of.write('virus\tquery\tquery_nt\tsubject\tsubject_leng\tEvalue\tLNVNRE\tIdentity\tAlign_query\tAlign_match\tAlign_subj\n')
		for index in cache[virus]:
			id=''
			of.write(virus+'\t')
			for i in xrange(11):
				if i ==0: continue
				out=linecache.getline(input, index+i).rstrip()
				if i==1:
					id='_'.join(linecache.getline(input, index+i).rstrip().split()[1:])
					of.write(id+'\t')
					if id in ids: continue
					print >>of2, '>', id+virus+'___'+base
				elif i==2:
					of.write(out.split()[1]+'\t')
					if id in ids: continue
					print >>of2, linecache.getline(input, index+i).strip().split()[1]
					ids.add(id)
				elif i==3:
					of.write(' '.join(out.split()[1:])+'\t')
				elif i==4 or i==5 or i==6 or i==7:
					of.write(out.split()[-1]+'\t')
				elif i==8 or i==9 or i==10:
					# of.write(out.split(' ')[-1]+'\t')
				# elif i==9:
					of.write('\n'+out)
			of.write('\n')
		of.close()
		of2.close()

if __name__ == '__main__': 
	input=sys.argv[1] #input of NR filtered blast output
	try: output=sys.argv[2] #sorted inputfile
	except: pass
	try: outfa=sys.argv[3]; 
	except: pass #sorted inputfileS
	cache, VE, NE, tax =CacheLines(input)
	OutputVirus(cache, input)
	printVirusSummary(cache, VE, NE, tax, input)