#!/usr/bin/env python

from collections import defaultdict
import operator
import sys
import os
import os.path
from os import listdir
from os.path import isfile, join

if __name__ == '__main__': 
	directory=sys.argv[1]
	wd=sys.argv[2]
	gooddb1 = [ f[0:(len(f)-4)] for f in listdir(wd) if (isfile(join(wd,f)) and f.endswith('.nal')) ] #multivolume db
	gooddb2 = [ f[0:(len(f)-4)] for f in listdir(wd) if (isfile(join(wd,f)) and f.endswith('_blastdb.nhr')) ] #singleVolume db
	gooddb = gooddb1+gooddb2
	os.chdir(wd)
	cmd ='blastdb_aliastool -dblist '+'\"'+' '.join(gooddb)+'\" -dbtype nucl -out '+directory+'_blastdb -title \"'+directory+'_blastdb\"'
	print cmd
	os.system(cmd)