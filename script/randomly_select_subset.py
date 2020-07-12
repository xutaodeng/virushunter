#Randomly select a subset of sequences (n, number sequnces) from fastq file

from Bio import SeqIO
import random

def subset_seq(filename):
    handle=open(filename)
    records=list(SeqIO.parse(handle,'fastq'))
    subset=random.sample(records, 5)
    SeqIO.write(subset, 'subset2.fastq','fastq')
    handle.close

subset_seq('good_quality.fastq')
    
    

