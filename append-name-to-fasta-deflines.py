#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
import sys

## This script will work after you run the command below from the "MERGED" directory that contains merged profile directories with the 
## SUMMARY-quality data structure that contains the MAG fasta files.  
outfile = open(sys.argv[1]+'_contigs.fa', 'w')
for seq in SeqIO.parse(open(sys.argv[1]+'temp.fa', 'r'), "fasta"):
    print(str(sys.argv[1])+"_"+str(seq.id))
    outfile.write(">"+str(sys.argv[1])+"_"+str(seq.id)+'\n'+str(seq.seq)+'\n')
