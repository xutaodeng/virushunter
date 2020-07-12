#!/usr/bin/env python

from collections import defaultdict
import operator
import sys
import os
import os.path
import linecache
import re, string

def plotpie(countfile):
	Rfile=countfile+'.R'
	of3 = open(Rfile, 'w')
	pieout=countfile+'.png'
	
	Rscript = "x=as.factor(c('Viridiplantae', 'NA', 'Fungi', 'Metazoa','Viruses', 'Phage', 'Archaea', 'Bacteria'))\n" +\
			"y=as.factor(  c('lightgreen',  'white','gray',  'blue',   'yellow',  'orange','pink'   ,'brown'))\n" +\
			' a<-read.table("'+countfile+'", header=F,  stringsAsFactors=F)\n'+ \
			'family <- as.numeric(a[[2]])\n'+\
			'lab<-a[[1]]\n'+\
			"lab[is.na(lab)]<-'NA'\n" +\
			'png("'+pieout+'", res=150, width=1000, height=1000)\n' +\
			'pie(family , col=as.character(droplevels(y[match(lab, x)])), labels=lab)\n'+\
			'dev.off()\n'
	of3.write(Rscript)
	of3.close()
	cmd='R CMD BATCH --quiet --vanilla '+Rfile
	print cmd
	os.system(cmd)

def outputHeader(of):
	head ='''
<html>
<head>
<link rel="stylesheet" type="text/css" href="../../tablestyle.css">
</head>
<script src="../../sorttable.js">
</script>
<script src="../../ajax_select.js">
</script>

'''
	of.write(head+'\n')

def printClark(countfile, htmlfile):
	of = open(htmlfile, 'w')
	f = open(countfile, 'r')

	outputHeader(of)
	of.write('<table class="sortable">\n')
	of.write('<tr><th>Order</th><th>Class</th><th>Family</th><th class="sorttable_numeric">HitCount</th><th>fasta</th></tr>\n')

	for line in f:
		if line.startswith('category,'): continue
		parts=line.strip().split(',')
		key='_'.join(parts[0:3]).replace('/', '_')
		key2=key.replace(',', '_')
		#print line.strip()
		try: order, clas, fam, count= parts
		except: 
			print parts
			sys.exit()
		fafile='fasta/'+os.path.basename(countfile)+'.'+key2+'.fa'
		of.write('<tr><td>'+order+'</td><td>'+clas+'</td><td>'+fam+'</td><td>'+count+'</td>')
		of.write('<td><a href=\"'+fafile+'\">fasta</a></td></tr>\n')
	of.write('</tbody></table></html>\n')
	of.close()
	f.close()

if __name__ == '__main__': 
	countfile1=sys.argv[1] #input of NR filtered blast output
	countfile=sys.argv[2] #input of NR filtered blast output
	htmlfile=sys.argv[3]
	printClark(countfile, htmlfile)
	#plotpie(countfile1)
