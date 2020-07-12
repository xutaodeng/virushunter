#!/usr/bin/env python
from collections import defaultdict
import os.path
import os
import string
import sys

all_virus=0
stats=defaultdict()
keys=['total_reads', 'raw_read_length', 'num_dup_reads', 'percent_dup', 'uni_reads', \
'clean_read_length', 'ave_bp_Overlap', 'nreads_Overlap', 'num_polyA_reads', 'percentHuman', 'percentBac', \
'num_tail_removed', 'tail_removed_average', '3prime_adaptors', '5prime_adaptors',\
'num_contig', 'num_contig_grater_60', 'num_seqs_to_blast', \
# 'Metazoa', \
# 'unmap',\
# 'Viridiplantae', \
# 'other_euk', \
# 'Fungi', \
# 'filter', \
# 'Viruses', \
# 'unclass', \
# 'Archaea', \
# 'Bacteria',\
# 'percent_virus_unmap_virus_unmap', \
'Sig_Vreads','n_hits', 'num_virus', 'n_mys',  'v_fam']

descriptions =["barcode: barcode", \
"total_reads: total number of reads for the barcode", \
"raw_read_length: raw read length",\
"num_dup_reads: number of duplicate clonal reads to be removed",\
"percentHuman: percentage of human reads", \
"percentBac: percentage of bacteria reads", \
# "num_polyA_reads: number of reads containing AAAAA",\
# "polyC: number of reads starting or ending with CCCCC",\
# "polyG: number of reads starting or ending with GGGGG",\
"percent_dup: percentage of duplicate reads to be removed",\
"uni_reads: number of unique reads", \
"clean_read_length: average cleaned read length",\
"ave_bp_Overlap: average overlap between the pair ends", \
"nreads_Overlap: number of paired reads showing overlap", \
"num_adaptors: number of adaptors",\
"num_tail_removed: number of tail removed by quality",\
#"tail_removed_average: Average bp trimed by quality", \
"num_contig: number of contigs formed", \
"num_contig_grater_100: number of contigs longer or equal to 100 bp", \
"num_seqs_to_blast: number of contigs or reads to blast virus", \
# "Metazoa: number of reads mapped to animals (near perfect alignment nucleotide)",\
# "unmap: unmapped to any species",\
# "Viridiplantae: mapped to plants",\
# "other_euk: mapped to other Euk than animals", \
# "Fungi: mapped to Fungi", \
# "filter: duplicates removed", \
# "Viruses: mapped to virus", \
# "unclass: mapped to unclassified species", \
# "Archaea: mapped to Archaea", \
# "Bacteria: mapped to Bacteria", \
# "percent_candidate_virus: percentage of candidate virus reads(virus, unmap, unclass)", \
"Sig_Vreads: number of sequences hit virus", \
"n_hits: number of virus hits passed NR filter", \
"num_virus: number of virus species (subspecies)", \
"n_mys: number of mysterious large contigs that doesn't hit to virus",\
#"vfam: virus hmmer alignments on the mysterious contigs",\
"v_fam: virus hmmer hits on the mysterious contigs"]

#sk= set(['raw_read_length', 'clean_read_length', 'ave_bp_Overlap', 'percentHuman', 'percentBac','percent_dup', 'percent_candidate_virus', 'tail_removed_average'])
sk= set(['raw_read_length', 'clean_read_length', 'percentHuman', 'percentBac','percent_dup'])

from virus_hunter import readSeeds2


def readstat(statfile='stats.logg'):
	global stats, all_virus
	f=open(statfile, 'r')
	for line in f:
		if '=' not in line:
			continue
		parts=line.split()
		try: id, name, value = parts[0], parts[-3], parts[-1]
		except: continue
		if 'all_blast_filter' in line: 
			all_virus= float(value)
		for key in stats.keys():
			if ('/'+key+'_' in id) or (key == id):
				try: stats[key][name].append(float(value))
				except: print  value, line
				break
	f.close()
	for key in stats.keys():
		for name in stats[key].keys():
			sum=0
			for value in stats[key][name]:
				sum+=value
			# if name!= 'percentHuman' and name!='ave_bp_Overlap' and (name in sk):
				# stats[key][name]=sum/2
			if name in sk:
				stats[key][name]=sum/len(stats[key][name])
			else:
				stats[key][name]=sum

def addallFile():
	global stats, all_virus
	out=defaultdict(float)
	for key in stats.keys():
		for name in stats[key].keys():
			out[name]+=stats[key][name]
	nlib=len(stats)
	for name in out.keys():
		if name in sk:
			out[name]/=nlib
	stats['all']=defaultdict()
	for name in out.keys():
		stats['all'][name]=out[name]
	stats['all']['num_virus'] = all_virus
	
def firstpage():
	f5 = open(wd+'/index.html', 'w')
	head ='''
	<html>
	<head>
		<style type="text/css" title="currentStyle">
			@import "../DataTables-1.9.4/media/css/demo_page.css";
			@import "../DataTables-1.9.4/media/css/demo_table.css";
		</style>
		<script type="text/javascript" language="javascript" src="../DataTables-1.9.4/media/js/jquery.js"></script>
		<script type="text/javascript" language="javascript" src="../DataTables-1.9.4/media/js/jquery.dataTables.js"></script>
		<script type="text/javascript" charset="utf-8">
			$(document).ready( function() {
			$('#example').dataTable( {
			"iDisplayLength": 1000,
			"sPaginationType": "full_numbers"
		} );
} )
		</script>
	</head>
	<table cellpadding="0" cellspacing="0" border="0" class="display" id="example">'''
	s=[]
	s.append('<thead><tr>')
	s.append('<th>barcode</th><th>Category table</th><th>NT Table</th><th>NT</th><th>pie table</th><th>pie figure</th><th>raw1</th><th>raw2</th><th>trimmed1</th><th>trimmed2</th>')
	for k in keys:
		s.append('<th>')
		s.append(k)
		s.append('</th>')
	s.append('</tr></thead>\n')
	head+=''.join(s)
	f5.write(head)
	f5.write('<tbody>')
	libs=stats.keys()
	libs = list(set(libs)-set(['all']))
	libs.append('all')
	for key in libs:
		outtxt = wd+'/blast_filter_out/'+key+'_blast_filter.txt'
		htmlFile =os.path.basename(os.path.splitext(outtxt)[0]+'.html')
		pietable ='pie/'+os.path.basename(os.path.splitext(outtxt)[0]+'.txt')
		mysfile ='mystery/'+key+'_m_complex_sort.fasta'
		mysfileout = 'mystery/'+key+'_m.fasta.out'
		mysfileout2 = 'mystery/'+key+'_m.fasta.out2'
		piepng ='pie/'+os.path.basename(os.path.splitext(outtxt)[0]+'.txt.png')
		clarktable ='clark/'+key+'.count'
		clarktable2 ='clark/'+key+'.html'
		clarkpng ='clark/'+key+'.count.png'
		#Pool21-25_raw_hist.png
		hist_raw_png1 ='hist/'+key+'_1_raw_hist.png'
		hist_raw_png2 ='hist/'+key+'_2_raw_hist.png'
		hist_clean_png1 ='hist/'+key+'_1_clean_hist.png'
		hist_clean_png2 ='hist/'+key+'_2_clean_hist.png'
		f5.write('<tr><td><a href="'+htmlFile+'">'+key+'</a></td>')
		
		f5.write('<td><a href="'+clarktable+'">'+key+'</a></td>')
		f5.write('<td><a href="'+clarktable2+'">'+key+'</a></td>')
		f5.write('<td><a href="'+clarkpng+'">'+'<img src="'+clarkpng+'" width="70" height="70"></a></td>')
		
		f5.write('<td><a href="'+pietable+'">'+key+'</a></td>')
		f5.write('<td><a href="'+piepng+'">'+'<img src="'+piepng+'" width="70" height="70"></a></td>')
		f5.write('<td><a href="'+hist_raw_png1+'">'+'<img src="'+hist_raw_png1+'" width="70" height="70"></a></td>')
		f5.write('<td><a href="'+hist_raw_png2+'">'+'<img src="'+hist_raw_png2+'" width="70" height="70"></a></td>')
		f5.write('<td><a href="'+hist_clean_png1+'">'+'<img src="'+hist_clean_png1+'" width="70" height="70"></a></td>')
		f5.write('<td><a href="'+hist_clean_png2+'">'+'<img src="'+hist_clean_png2+'" width="70" height="70"></a></td>')
		
		for k in keys:
			out='<td>'
			if k=='n_mys': out+='<a href="'+mysfile+'">'
			#if k=='vfam': out+='<a href="'+mysfileout+'">'
			if k=='v_fam': out+='<a href="'+mysfileout2+'">'
			try: out+=str(int(stats[key][k]))
			except: out+='-'
			if k=='n_mys': out+='</a>'
			#if k=='vfam': out+='</a>'
			if k=='vfam2': out+='</a>'
			out+='</td>'
			f5.write(out)
		f5.write('</tr>\n')
	
	# outtxt = wd+'/blast_filter_out/'+'all'+'_blast_filter.txt'
	# htmlFile =os.path.basename(os.path.splitext(outtxt)[0]+'.html') #last file
	# f5.write('<a href="'+htmlFile+'">'+htmlFile+'</a><br>')
	f5.write('</tbody>')
	f5.write('</table>')
	for des in descriptions:
		f5.write("<br>"+des)
	f5.write('</html>')
	f5.close()

if __name__ == "__main__":
	reAssemble = sys.argv[1]
	stats, seeds=readSeeds2(reAssemble)
	wd = os.path.abspath(os.path.dirname('.'))
	print wd
	readstat()
	addallFile()
	firstpage()