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
<head>
<link rel="stylesheet" type="text/css" href="../tablestyle.css">
</head>
<script src="../sorttable.js">
</script>
<script src="../ajax_select.js">
</script>
E-value output range: from <input type="text" id="evalue1" name="evalue1" value="0"> to <input type="text" id="evalue2" name="evalue2" value="1"><br>
<button type="button" onclick="update_Evalue_range()">Update Table </button><br>
<br>
<input type="checkbox" name='checkall' value='checkall' onclick="checkAll(this)">
Check All <button onclick="exportData()" >Export Selected</button> <br>

 <div id="txtHint"> </div>

'''
#<tr><th class="sorttable_nosort"><th>Virus</th><th class="sorttable_numeric">Low_V_Evalue</th><th class="sorttable_numeric">Low_NR_Evalue</th><th class="sorttable_numeric">Hits</th><th>fasta</th>


	of.write(head+'\n')

# def hitHeader(of):
	# head ='''
# <html>
# <head>
# <link rel="stylesheet" type="text/css" href="../../tablestyle.css">
# </head>
# <script src="../../sorttable.js">
# </script>
# <table class="sortable">
# <tr><th class="sorttable_nosort">Virus</th><th>Query</th><th>Pair Hit</th><th>Evalue</th><th>LNVNRE</th><th>Identity</th><th class="sorttable_nosort">Alignment</th><th>Query_nt</th></tr>
# '''
	# of.write(head+'\n')
	
def hitHeader(of):
	head ='''
	<html>
	<head>
		<style type="text/css" title="currentStyle">
			@import "../../DataTables-1.9.4/media/css/demo_page.css";
			@import "../../DataTables-1.9.4/media/css/demo_table.css";
		</style>
		<script type="text/javascript" language="javascript" src="../../DataTables-1.9.4/media/js/jquery.js"></script>
		<script type="text/javascript" language="javascript" src="../../DataTables-1.9.4/media/js/jquery.dataTables.js"></script>
		<script type="text/javascript" charset="utf-8">
			$(document).ready( function() {
			$('#example').dataTable( {
			"iDisplayLength": 500,
			"sPaginationType": "full_numbers"
		} );
} )
		</script>
	</head>
	<table cellpadding="0" cellspacing="0" border="0" class="display" id="example">
	<thead><tr><th>Virus</th><th>Query</th><th>Pair Hit</th><th>Evalue</th><th>LNVNRE</th><th>Identity</th><th>Alignment</th><th>Query_nt</th><th>Pair_nt</th></tr></thead><tbody>
'''
	
	# head ='''
# <html>
# <head>
# <link rel="stylesheet" type="text/css" href="../../tablestyle.css">
# </head>
# <script src="../../sorttable.js">
# </script>
# <table class="sortable">
# <tr><th class="sorttable_nosort">Virus</th><th>Query</th><th>Pair Hit</th><th>Evalue</th><th>LNVNRE</th><th>Identity</th><th class="sorttable_nosort">Alignment</th><th>Query_nt</th></tr>
# '''
	of.write(head+'\n')

def CacheFASTA(fname, hit_pairs): #for retrieving pair sequences
	cachecombine={}
	f = open(fname, 'r')
	i=0
	start, end = 0,0
	header= None
	print 'caching fasta'
	for line in f:
		i+=1
		if i%10000000==0: print 'caching fasta line', i
		if line.strip().startswith('>'):
			end=i-1
			if header!=None and (header in hit_pairs): 
				cachecombine[header] = (start, end)
			line=line[1:] #remove >
			if line.startswith('@'): #read
				parts = line.strip().split('_')
				header= parts[0]+'_'+parts[1]
			else: header='' #contig, ignore
			start=i+1
	if header!=None and (header in hit_pairs): cachecombine[header] = (start, i)
	print 'cacheFASTA done'
	return cachecombine

def getPairSeq(header, combine, cachecombine):
	seq=[]
	try: start, end = cachecombine[header]
	except: return ''
	for i in xrange(start, end+1):
		seq.append(linecache.getline(combine, i).strip())
	return ''.join(seq)

def CacheLines(input): 
	#print 'cache line'
	cache=defaultdict(list)
	c2, c5, c10 = defaultdict(int), defaultdict(int), defaultdict(int)
	VE={}
	NE={}
	tax={}
	qlen={}
	hit_pairs=set([])
	f = open(input, 'r')
	i=0
	for line in f:
		#if i%10000000==0: print i
		i+=1
		if i%11 == 2 and line.startswith('@'):
			parts=line.strip().split('_')
			if parts[1]=='1': hit_pairs.add(parts[0]+'_'+'2') #e.g s12576_1
			elif parts[1]=='2': hit_pairs.add(parts[0]+'_'+'1')
		if i%11 == 3:
			try: query = line.strip().split()[1]
			except: 
				print i, line
				sys.exit(1)
		if i%11 == 4:
			taxonomy=line.strip().split()[-1]
			try: virus=taxonomy.split(':')[0].split('$')[1]
			except: print i, line; sys.exit()
			cache[virus].append(i-3)
			tax[virus]=taxonomy
			if not qlen.has_key(virus):
				qlen[virus]=len(query)
			elif qlen.has_key(virus) and qlen[virus]<len(query):
				qlen[virus]=len(query)
		elif i%11==6: #veval
			ve=float(line.split()[-1])
			if ve<=1E-2: 
				c2[virus]+=1
			if ve<=1E-5: 
				c5[virus]+=1
			if ve<=1E-10: 
				c10[virus]+=1
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
	#print 'cache line done'
	return cache, VE, NE, tax, qlen, c2, c5, c10, hit_pairs

def printVirusSummary(cache, VE, NE, tax, input, qlen, c2, c5, c10, vcounts, barcodes):
	try: os.mkdir(os.path.dirname(input)+'/pie/')
	except: pass
	try: os.mkdir(os.path.dirname(input)+'/table/')
	except: pass
	base=os.path.basename(input)
	of2 = open(os.path.dirname(input)+'/pie/'+base, 'w')
	of3 = open(os.path.dirname(input)+'/table/'+base, 'w')
	of = open(os.path.splitext(input)[0]+'.html', 'w')
	f4=os.path.splitext(input)[0]+'.xls'
	ff4 = os.path.basename(f4)

	keys=sorted(cache.keys())
	outputHeader(of)
	of.write('<a href=\"'+ff4+'\">excel_download</a>\n')
	of.write('<a href=\"hitTable_e-2.xls\">Table_Evalue10-2</a>\n')
	of.write('<a href=\"hitTable_e-5.xls\">Table_Evalue10-5</a>\n')
	of.write('<a href=\"hitTable_e-10.xls\">Table_Evalue10-10</a>\n')
	of.write('<table class="sortable">\n')
	of.write('<tr><th class="sorttable_nosort"></th><th>Category</th><th>Family</th><th>Genus</th><th>Species</th><th class="sorttable_numeric">Max_Contig</th><th class="sorttable_numeric">Low_V_Evalue</th><th class="sorttable_numeric">Low_NR_Evalue</th><th class="sorttable_numeric">Hits</th><th>fasta</th>')
	for barcode in barcodes:
		of.write('<th>'+barcode+'</th>')
	of.write('</tr>\n')

	# of4.write('Category\tFamily\tGenus\tSpecies\tMax_Contig\tLow_V_Evalue\tLow_NR_Evalue\tHits\tfasta')
	# for barcode in barcodes:
		# of4.write('\t'+barcode)
	# of4.write('\n')
	numvirus=0
	pie=defaultdict(int)
	for virus in keys:
		try: 
			species, genus, family, cat = tax[virus].split(':')[0:4]
			cat=cat.split('$')[1]
			family=family.split('$')[1]
			genus=genus.split('$')[1]
			species=species.split('$')[1]
		except: cat,species, genus, family = 'NA', 'NA', 'NA', 'NA'
		base=os.path.basename(input)
		# virname='_'.join(virus.strip('[').strip(']').replace('/','_').split())[0:40]
		# virname = re.sub('[\W_]+', '_', virname)
		#virname=species.replace('/','_').replace('*', '_')
		virname=''
		for e in species:
			if e.isalnum(): virname+=e
			else: virname+='_'
		label = base+'_'+virname
		f1 = 'aln/'+label+'.html'
		f2 = 'fasta/'+label+'.fa'
		of.write('<tr><td><input type="checkbox" value="'+label+'"> </td>')
		#of.write('<td>'+virname+'</td><td>'+str(VE[virus])+'</td><td>'+str(NE[virus])+'</td><td><a href=\"'+f1+'\">'+str(len(cache[virus]))+'</a></td><td>')
		of.write('<td>'+cat+'</td><td>'+family+'</td><td>'+genus+'</td><td>'+species+'</td><td>'+str(qlen[virus])+'</td><td>'+str(VE[virus])+'</td><td>'+str(NE[virus]))
		of.write('</td><td><div id=\"'+label+'_hit\"><a href=\"'+f1+'\">'+str(len(cache[virus]))+'</a></div></td><td>')
		of.write('<div id=\"'+label+'_fa\"><a href=\"'+f2+'\">fasta</a></div></td>')
		
		# of4.write(cat+'\t'+family+'\t'+genus+'\t'+species+'\t'+str(qlen[virus])+'\t'+str(VE[virus])+'\t'+str(NE[virus]))
		# of4.write('\t'+f1+'\t'+str(len(cache[virus]))+'\t'+f2+'\n')
		
		for barcode in barcodes:
			try: vc=vcounts[virname][barcode]
			except: vc=0
			of.write('<td>'+str(vc)+'</td>')
		of.write('</tr>\n')
		of3.write(cat+'\t'+family+'\t'+genus+'\t'+species+'\t'+str(c2[virus])+'\t'+str(c5[virus])+'\t'+str(c10[virus])+'\n')
		pie[family]+=len(cache[virus])
		numvirus+=1
	of.write('</tbody></table></html>\n')
	of.close()
	os.system('cp ' +os.path.splitext(input)[0]+'.html '+f4)
	for (family, count) in pie.items():
		of2.write(family+' '+str(count)+'\n')
	of2.close()
	of3.close()
	print input, 'num_virus =', numvirus


# def plotpie(input):
	# base=os.path.basename(input)
	# of3 = open(os.path.dirname(input)+'/pie/'+base+'.R', 'w')
	# piein=os.path.dirname(input)+'/pie/'+base
	# pieall=os.path.dirname(input)+'/pie/all_blast_filter.txt'
	# pieout=os.path.dirname(input)+'/pie/'+base+'.png'
	
	# Rscript =' aa<-read.table("'+piein+'", header=F,  stringsAsFactors=F)\n'+ \
			# ' s<-read.table("'+pieall+'", header=F,  stringsAsFactors=F)\n'+ \
			# ' a<-merge(s, aa, by="V1", all=T)\n'+\
			# ' a[[3]][is.na(a[[3]])]<-0\n'+ \
			# # 'b<-a[ a[[2]]/sum(a[[2]])>0.01,]\n' +\
			# # 'd<-a[ a[[2]]/sum(a[[2]])<=0.01,]\n' +\
			# # 'x<-rbind(b, c("Others", sum(d[[2]])))\n' +\
			# # 'family <- as.numeric(x[[2]])\n' +\
			# # 'lab<-x[[1]]\n' +\
			# 'family <- as.numeric(a[[3]])\n'+\
			# 'lab<-a[[1]]\n'+\
			# 'lab[family==0]<-NA\n'+\
			# 'png("'+pieout+'", res=150, width=1000, height=1000)\n' +\
			# 'pie(family, , col=rainbow(length(family)), labels=lab)\n'+\
			# 'dev.off()\n'
	# of3.write(Rscript)
	# of3.close()
	# cmd='R CMD BATCH --quiet --vanilla '+os.path.dirname(input)+'/pie/'+base+'.R'
	# print cmd
	# os.system(cmd)


def detectPair(cache, input):
	#print 'start pairing'
	keys=sorted(cache.keys())
	#print 'len( virus keys)', len(keys)
	allpairs={}
	vv=0
	for virus in keys:
		vv+=1
		#if vv%1000==0: print 'num_virus processed', vv
		allpairs[virus]=defaultdict(set)
		for index in cache[virus]:
			i=1
			id=linecache.getline(input, index+i).strip()
			if id.startswith('@'):
				position, pair = id.split('_')[0:2]
				allpairs[virus][position].add(pair)
	return allpairs

def OutputVirus(cache, input, allpairs, combine, cachecombine, base, cwd):
	keys=sorted(cache.keys())
	try: os.mkdir(os.path.dirname(input)+'/aln/')
	#print 'dir', os.path.dirname(input)+'/aln/'
	except: pass
	try: os.mkdir(os.path.dirname(input)+'/fasta/')
	except: pass
	countfile=os.path.dirname(input)+'/aln/viralcounts.txt'
	fc=open(countfile, 'w')
	fc.close()
	try: os.mkdir(os.path.dirname(input)+'/tmp/')
	except: pass
	counter=0
	for virus in keys:
		try: 
			species, genus, family, cat = tax[virus].split(':')[0:4]
			cat=cat.split('$')[1]
			family=family.split('$')[1]
			genus=genus.split('$')[1]
			species=species.split('$')[1]
		except: cat,species, genus, taxonomy = 'NA', 'NA', 'NA', 'NA'
		#ids=set()
		base=os.path.basename(input)
		#virname=species.replace('/','_')
		virname=''
		for e in species:
			if e.isalnum(): virname+=e
			else: virname+='_'
		label = base+'_'+virname
		of = open(os.path.dirname(input)+'/aln/'+label+'.html', 'w')
		of2 = open(os.path.dirname(input)+'/tmp/'+label+'.fa', 'w')
		#of.write('virus\tquery\tquery_nt\tsubject\tsubject_leng\tEvalue\tLNVNRE\tIdentity\tAlign_query\tAlign_match\tAlign_subj\n')
		hitHeader(of)
		for index in cache[virus]:
			id=''
			query_nt=''
			pair_nt=''
			outrow=['<tr><td>'+virus+'<td>']
			for i in xrange(11):
				if i ==0: 
					continue
				out=linecache.getline(input, index+i).rstrip()
				#seqid='@s'+str(lineno)+'_'+pair_end+'_'+label
				if i==1:
					id=out.strip()
					#print id
					if id.startswith('@'):
						outid=id
						paired=1
						try:
							position, pair = id.split('_')[0:2]
							if pair=='1': pp='2'
							elif pair=='2': pp='1'
							else: pass#print 'err pair', pair
							#key = id.split('_',2)[2]
							header=position+'_'+pp #+'_'+key
							pair_nt= getPairSeq(header, combine, cachecombine)
							paired=len(allpairs[virus][position])
						except: pass
					else:
						outid='contig_'+'_'.join(out.split())
						paired=0
					outrow.append(outid+'<td>')
					if paired==2: outrow.append('Y<td>')
					else: outrow.append('-<td>')
					#if id in ids: continue
					print >>of2, '>'+id+' '+virname
				elif i==2:
					query_nt =out.strip().split()[1]
					#if id in ids: continue
					print >>of2, linecache.getline(input, index+i).strip().split()[1]
					#ids.add(id)
				elif i==3 or i==4:
					pass #outrow.append(' '.join(out.split()[1:])+'<td>')
				elif i==5 or i==6 or i==7:
					outrow.append(out.strip().split()[-1]+'<td>')
				elif i==8 or i==9:
					outrow.append(out[11:]+'<BR>')
					# of.write(out.split(' ')[-1]+'\t')
				elif i==10:
					outrow.append(out[11:]+'<td>')
					nn=len(query_nt)/3
					qt = query_nt[0:nn]+'<BR>'+query_nt[nn:2*nn]+'<BR>'+query_nt[2*nn:]
					outrow.append(qt+'<td>')
					nn=len(pair_nt)/3
					pair_qt = pair_nt[0:nn]+'<BR>'+pair_nt[nn:2*nn]+'<BR>'+pair_nt[2*nn:]
					outrow.append(pair_qt+'</tr>')
					of.write(''.join(outrow)+'\n')
			#of.write()
		of.write('</table></html>\n')
		of.close()
		of2.close()
		fain=os.path.dirname(input)+'/tmp/'+label+'.fa'
		faout=os.path.dirname(input)+'/fasta/'+label+'.fa'
		os.system(dirscr+'faSort.py '+fain+' '+faout+' False') #sort fa by length
		counter+=1
		if combine =='--': #'all' barcode
			#if counter>5: break
			#base is all_blast_filter.txt
			cmd = dirscr+'viralCount.py '+faout+' '+countfile+' '+cwd+' '+virname
			#print cmd
			os.system(cmd) #sort fa by length
	vcounts={}
	barcodes=set([])
	if combine =='--': #'all' barcode
		vcounts, barcodes=readCount(countfile)
		barcodes=list(barcodes)
		barcodes=sorted(barcodes)
	return vcounts, barcodes

def readCount(countfile):
	vcounts={}
	barcodes=set([])
	f=open(countfile, 'r')
	for line in f:
		parts=line.strip().split()
		barcode, virname, count=parts
		if not vcounts.has_key(virname):
			vcounts[virname]={}
		vcounts[virname][barcode]=count
		barcodes.add(barcode)
	f.close()
	return vcounts, barcodes

if __name__ == '__main__': 
	input=sys.argv[1] #input of NR filtered blast output
	cache, VE, NE, tax, qlen, c2, c5, c10,hit_pairs =CacheLines(input)
	
	#plotpie(input)
	combine=sys.argv[2] #_r file combined for extracting , no len filter fasta
	dirscr=sys.argv[3]
	try: cachecombine=CacheFASTA(combine, hit_pairs) #combined fasta file cache
	except: cachecombine ={}
	base=sys.argv[4] #sorted inputfile
	cwd=sys.argv[4] #cwd
	try: outfa=sys.argv[6]; 
	except: pass #sorted inputfileS
	allpairs={}
	#allpairs = detectPair(cache, input)
	vcounts, barcodes=OutputVirus(cache, input, allpairs, combine, cachecombine, base, cwd) #individual read page
	printVirusSummary(cache, VE, NE, tax, input, qlen, c2, c5, c10, vcounts, barcodes) #summary page
