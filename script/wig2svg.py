#!/usr/bin/env python
import sys
import os
from collections import defaultdict
from operator import itemgetter, attrgetter
import math

wigfile = sys.argv[1]
svgfile=sys.argv[2]


def readWig4(wigfile):
	f=open(wigfile, 'r')
	wigs=[]
	for line in f:
		parts=line.strip().split()
		ref, start, end, wig=parts
		wigs.append((int(start), int(end), float(wig)))
	f.close()
	return wigs
	
def readWig2(wigfile):
	f=open(wigfile, 'r')
	wigs=[]
	prevPos, prevValue=0, 0
	header=f.readline()
	for line in f:
		parts=line.strip().split()
		start, wig=parts
		if prevPos> 0: 
			wigs.append((prevPos, int(start), float(prevValue)))
		prevPos, prevValue=int(start), float(wig)
		#wigs.append((int(start), int(end), float(wig)))
	wigs.append((prevPos, prevPos+1, float(prevValue)))
	f.close()
	return wigs

def print_SVG(wigs, svgfile):
	header ='''
	<!DOCTYPE html>
	<html>
	<body>
	'''
	of=open(svgfile, 'w')
	of.write('<br>\n')
	label=svgfile.rsplit('.',1)[0].rsplit('/', 1)[-1]
	scale = 10
	yscale=20
	body=['<svg height="200" width="1200">']
	body.append('<text x="200" y="160" fill="green">'+label+'</text>')
	body.append('<line x1="1000" y1="0" x2="1000" y2="120" style="stroke:rgb(255,0,0);stroke-width:3" />')
	body.append('<line x1="0" y1="130" x2="1000" y2="130" style="stroke:rgb(255,0,0);stroke-width:3" />')
	for depth in [0, 10, 100, 1000, 10000, 100000, 1000000]:
	#for i in xrange(6):
	#	depth=math.pow(10, i)-1
		ypos=str(math.log(depth+1, 10)*yscale)
		ypos2=str(math.log(depth+1, 10)*yscale+10)
		body.append('<text x="1010" y="'+ypos2+'" fill="red">'+str(int(depth))+'</text>')
		body.append('<line x1="1000" y1="'+ypos+'" x2="1010" y2="'+ypos+'" style="stroke:rgb(255,0,0);stroke-width:3" />')
	for i in xrange(10):
		xpos = i*1000/scale
		x=i*1000
		body.append('<line x1="'+str(xpos)+'" y1="130" x2="'+str(xpos)+'" y2="120" style="stroke:rgb(255,0,0);stroke-width:3" />')
		body.append('<text x="'+str(xpos)+'" y="145" fill="red">'+str(i)+'k'+'</text>')
	for (start, end, wig) in wigs:
		x=str(float(start)/scale)
		width=str(float(end-start)/scale)
		height=str(math.log(wig+1, 10)*yscale)
		body.append('<rect x="'+x+'" y="0" width="'+width+'" height="'+height+'" style="fill:rgb(0,0,255);stroke-width:1;stroke:rgb(0,0,255)" />')
	

	# trailer ='''</svg>
	# </body>
	# </html>
	# '''
	body.append('</svg>')
	print >> of, '\n'.join(body)

#wigs=readWig4(wigfile) #4 column format
wigs=readWig2(wigfile) #2 column format
print_SVG( wigs, svgfile)