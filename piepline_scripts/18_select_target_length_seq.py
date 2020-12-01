#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
author: jin
date: April.13,2018
function: this code is used to rename the names of every seqences and delete empty sequences in fasta file.
usage: python ~ <fasta file> <Output file> <target length> 
"""
import sys, os
from Bio import SeqIO
if len(sys.argv) != 4:
	print "Usage: python ~ <fasta file> <Output file> <target length>"
	sys.exit()


outfile = open(sys.argv[2], 'w')
###########################################################################
#set counters
c1=0
c2=0
for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
    c1 += 1
    if len(seq_record) > int(sys.argv[3]):
        outfile.write(">" + seq_record.id.strip() + "\n" + str(seq_record.seq) + "\n")
        c2 += 1
print "The total number of sequences is %s, and there are %s greater than %s." %(c1,c2,sys.argv[3]) 

outfile.close()
