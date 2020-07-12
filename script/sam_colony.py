#!/usr/bin/env python
import operator
import sys
#print 'running...'
f=open(sys.argv[1],'r')
out=open(sys.argv[2],'w')
out2=open(sys.argv[3],'w')
out3=open(sys.argv[4],'w')

b=0 # clonal
d={}
r=0 # aligned
u=0 # unaligned

for line in f:
    if line[0]=='@':  # sam file header line
        print>>out,line.strip()
        continue
    current=line.strip().split()
    if line.strip().split()[2]=='*': # unaligned
        u+=1
        continue
    r+=1
    if r%10000000==0:
        print r
    
    flag=bin(int(current[1]))
    
    try: strand=flag[-5]
    except: 
        #print 'flag', flag
        strand = '0'
    kee = current[2]+'_'+current[3]+'_'+strand
    if d.has_key(kee):
        b+=1
        #print>>out2,'\t'.join(current)
        d[kee]+=1
    else:
        d[kee]=0
        print>>out,'\t'.join(current)

print >>out2, 'unaligned reads in sam: ', u
print >>out2, 'total aligned reads in sam: ', r
print >>out2, "total unique reads: ", len(d)
print >>out2,'colony reads in sam: ', b
print >>out2,'colony reads percentage:', str(b/float(r)*100)+'%'

for key in d.keys():
    if d[key]<10:
        del d[key]
print >>out2, 'total unique reads removing low freq (<10) alignments', len(d)

sorted_d=sorted(d.iteritems(),key=operator.itemgetter(1), reverse=True)
d={}

for i in xrange(len(sorted_d)):
    print >>out3,sorted_d[i][0], sorted_d[i][1]
