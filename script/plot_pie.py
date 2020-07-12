#!/usr/bin/env python

from collections import defaultdict
import operator
import sys
import os
import os.path
import linecache
import re, string

def plotpie(input):
	base=os.path.basename(input)
	of3 = open(os.path.dirname(input)+'/'+base+'.R', 'w')
	piein=os.path.dirname(input)+'/'+base
	pieall=os.path.dirname(input)+'/all_blast_filter.txt'
	pieout=os.path.dirname(input)+'/'+base+'.png'
	
	Rscript =' aa<-read.table("'+piein+'", header=F,  stringsAsFactors=F)\n'+ \
			' s<-read.table("'+pieall+'", header=F,  stringsAsFactors=F)\n'+ \
			' a<-merge(s, aa, by="V1", all=T)\n'+\
			' a[[3]][is.na(a[[3]])]<-0\n'+ \
			'family <- as.numeric(a[[3]])\n'+\
			'lab<-a[[1]]\n'+\
			'lab[family==0]<-NA\n'+\
			'png("'+pieout+'", res=150, width=1000, height=1000)\n' +\
			'pie(family, , col=rainbow(length(family)), labels=lab)\n'+\
			'dev.off()\n'
	of3.write(Rscript)
	of3.close()
	cmd='R CMD BATCH --quiet --vanilla '+os.path.dirname(input)+'/'+base+'.R'
	print cmd
	os.system(cmd)

if __name__ == '__main__': 
	input=sys.argv[1] #input of NR filtered blast output
	plotpie(input)
