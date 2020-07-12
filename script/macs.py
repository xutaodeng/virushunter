#!/usr/bin/env python
import os

#pvalue='1e-5' #defaut for macs
pvalue='1e-9' #defaut for super enhancer

f=open('input.sams', 'r')
of=open('macs.sh','w')
of2=open('enhancer.sh','w')
of3=open('peak_readCounts.sh','w')

servers=['bsidna3', 'bsidna4', 'bsidna5','bsidna6','bsidna7','bsidna8','bsidna9','bsidna10','bsidna11',\
'bsidna12','bsidna13','bsidna14','bsidna15','bsidna16','bsidna17','bsidna18','bsidna19', 'bsidna20', \
'bsidna21', 'bsidna22', 'bsidna23', 'bsidna24','bsidna25', 'bsidna26', 'bsidna27','bsidna28','bsidna29', 'bsidna30',\
'bsidna32', 'bsidna33', 'bsidna34', 'bsidna35', 'bsidna36']

job=0
wd=os.getcwd()
print 'wd', wd
for line in f:
	id=line.strip().split('/')[-1].split('.')[0].split('_',1)[1]
#	of2.write('ssh '+servers[job%len(servers)]+ ' mv bowtie* '+wd+'\n')
	if job%2==0:
		treat=line.strip()#[1:]
	else:
		control=line.strip()#[1:]
		treat, control = control, treat
		of.write('ssh '+servers[job%len(servers)]+ ' " cd '+ wd+ \
		' && /mnt/cluster/tools/MACS/bin/macs14 -p '+pvalue+' --format SAM -t '+ \
		treat+' -c '+control+' -n '+ id+' " &\n')
		of2.write('ssh '+servers[job%len(servers)]+ ' " cd '+ wd+ \
		' && superEnhancer.py '+id+'_peaks.bed '+ id+'.enhancer " &\n')
		of3.write('ssh '+servers[job%len(servers)]+ ' " cd '+ wd+ \
		' && peak_readCounts.py '+id+'.enhancer '+treat+' '+ id+'.treat.count '+control+' " &\n')
		# of3.write('ssh '+servers[job%len(servers)]+ ' " cd '+ wd+ \
		# ' && peak_readCounts.py '+id+'.enhancer '+control+' '+ id+'.control.count " &\n')
		
	job+=1
#this is for the last sample without input control, treat only
# of.write('ssh '+servers[job%len(servers)]+ ' " cd '+ wd+ \
		# ' && /mnt/cluster/tools/MACS/bin/macs14 --format SAM -t '+treat+' -n '+ id+' " &\n')
# of2.write('ssh '+servers[job%len(servers)]+ ' " cd '+ wd+ \
		# ' && superEnhancer.py '+id+'_peaks.bed '+ id+'.enhancer " &\n')
# of3.write('ssh '+servers[job%len(servers)]+ ' " cd '+ wd+ \
		# ' && peak_readCounts.py '+id+'.enhancer '+treat+' '+ id+'.treat.count " &\n')

f.close()
of.close()
of2.close()
of3.close()

of4=open('enhancer_run.sh', 'w')
of4.write('source macs.sh\nwait\n')
of4.write('source enhancer.sh\nwait\n')
of4.write('source peak_readCounts.sh\nwait\n')
of4.close()

