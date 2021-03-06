#!/usr/bin/env python
from collections import defaultdict
import sys
import os
import os.path
import itertools
import gzip
import linecache
import gc
from operator import itemgetter
#import sqlite3
import gzip

count=0
total=0

# def RestoreOrderDB(tmpfq, outfq):
	# try: os.remove(tmpfq+'.db')
	# except: pass
	# #conn = sqlite3.connect(tmpfq+'.db')
	# conn = sqlite3.connect(":memory:")
	# c = conn.cursor()

	# # Create table
	# c.execute('''CREATE TABLE ordertable 
               # (lineNo integer, ind integer)''')

	# # c.execute('''CREATE TABLE stocks
             # # (date text, trans text, symbol text, qty real, price real)''')
			 
	# f=open(tmpfq, 'r')
	# i=0
	# for line in f:
		# if i%1000000==0: print i
		# i+=1
		# lineNo = int(line.strip().split()[0])
		# c.execute("INSERT INTO ordertable VALUES ("+str(lineNo)+","+str(i)+")")
	# f.close()
	# conn.commit()

	# of=open(outfq, 'w')
	# j=0
	# for row in c.execute('SELECT * FROM ordertable ORDER BY lineNo'):
		# lineNo, i = row[0], row[1]
		# if j%1000000==0: print 'j',j
		# j+=1
		# line = linecache.getline(tmpfq, i)
		# Number, content = line.strip().split()
		# of.write(content+'\n')
	# of.close()

	# conn.close()

def RestoreOrder(tmpfq, outfq):
	f=open(tmpfq, 'r')
	order = []
	i=0
	for line in f:
		if i%1000000==0: print i
		#if i>10000000: break
		lineNo = int(line.strip().split()[0])
		#print lineNo, line
		i+=1
		if i%4==1: order.append((lineNo,i))
	f.close()
	order.sort(key=itemgetter(0))
	of=open(outfq, 'w')
	j=0
	for (lineNo, i) in order:
		try:
			#if j%1000000==0: print j
			j+=1
			line = linecache.getline(tmpfq, i)
			Number, content = line.strip().split()
			of.write(content+'\n')
			line = linecache.getline(tmpfq, i+1)
			Number, content = line.strip().split()
			of.write(content+'\n')
			line = linecache.getline(tmpfq, i+2)
			Number, content = line.strip().split()
			of.write(content+'\n')
			line = linecache.getline(tmpfq, i+3)
			Number, content = line.strip().split()
			of.write(content+'\n')
		except:
			pass#print 'exception', lineNo, i, j, line
	of.close()

# def RestoreOrder(tmpfq, outfq):
	# f=open(tmpfq, 'r')
	# order = []
	# i=0
	# for line in f:
		# if i%1000000==0: print i
		# i+=1
		# lineNo = int(line.strip().split()[0])
		# order.append((lineNo,i))
	# f.close()
	# order.sort(key=itemgetter(0))
	# of=open(outfq, 'w')
	# for (lineNo, i) in order:
		# line = linecache.getline(tmpfq, i)
		# Number, content = line.strip().split()
		# of.write(content+'\n')
	# of.close()


def removedup(inputfq, outputfq):
	#print inputfq
	global count, total
	#print 'processing..', inputfq
	if inputfq.endswith('.gz'):
		f=gzip.sopen(inputfq, 'rb')
	else:
		f=open(inputfq, 'r')
	of=open(outputfq, 'a')
	i=0
	keyset=set()
	for line in f:
		i+=1
		if i%4==1:
			id = line.strip()
		elif i%4==2:
			readNo, read = line.strip().split()
		elif i%4==3:
			id2 = line.strip()
		elif i%4==0:
			total+=1
			qualNo, qual = line.strip().split()
			if read[0:50] in keyset: read='A'; qual='A'; count+=1
			print >>of, '\n'.join([id,readNo+' '+read,id2,qualNo+' '+qual])
			keyset.add(read[0:50])
	f.close()
	of.close()

def split(filename):
	if filename.endswith('.gz'):
		f=gzip.open(filename, 'r')
	else:
		f=open(filename, 'r')
	comb=itertools.product(['A','C','G','T'], repeat=4)
	keydict={'unknown':open(filename+'_unknown','w')}
	for it in comb: 
		key=''.join(it)
		keydict[key]=open(filename+'_'+key,'w')

	i=0
	for line in f:
		i+=1
		if i%4==1:
			id = line.strip().split()[0]
		elif i%4==2:
			read = line.strip()
			readstart=read[0:4]
		elif i%4==3:
			id2 = line.strip()
		elif i%4==0:
			qual = line.strip()
			if keydict.has_key(readstart): print >>keydict[readstart], '\n'.join([str(i-3)+' '+id,str(i-2)+' '+read,str(i-1)+' '+id2,str(i)+' '+qual])
			else: print >>keydict['unknown'], '\n'.join([str(i-3)+' '+id,str(i-2)+' '+read,str(i-1)+' '+id2,str(i)+' '+qual])
	f.close()
	for key in keydict.keys():
		keydict[key].close()

def cleanTempFiles():
	comb=itertools.product(['A','C','G','T'], repeat=4)
	keydict={}
	for it in comb: 
		key=''.join(it)
		os.remove(filename+'_'+key)
	os.remove(filename+'_unknown')
	os.remove(filename+'.tmp')

if __name__ == '__main__':
	try:
		filename=sys.argv[1]
		outfile=sys.argv[2]
		tmpfq=filename+'.tmp'
		of=open(tmpfq, 'w') #clear this output file
		of.close()
		split(filename)
		comb=itertools.product(['A','C','G','T'], repeat=4)
		removedup(filename+'_unknown', tmpfq)
		for it in comb: 
			key=''.join(it)
			removedup(filename+'_'+key, tmpfq)
		gc.collect()

		RestoreOrder(tmpfq, outfile)
		cleanTempFiles()
	except (KeyboardInterrupt, SystemExit):
		cleanTempFiles()
	except:
		print 'usage: dedup.py input.fastq output.deduped.fastq'
		sys.exit()
	print filename, 'num_dup_reads =', count
	print filename, 'percent_dup =', round(count/float(total)*100, 2) 
	#print filename, 'total', total, 'dup', count, 'percentDup', round(count/float(total)*100, 2) 
