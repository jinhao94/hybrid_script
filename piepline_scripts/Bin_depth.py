#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
author: jin
date: 2019.1.29
function: statistics depth of each bin 
"""

import sys, os
import argparse

parser = argparse.ArgumentParser(description='Calculate the normalized depth of each bin')


parser.add_argument("-b", '--bin', required=True, help='Bin group file (from Fasta_to_Scaffolds2Bin.sh)')
parser.add_argument("-d", '--depth' , default= True , help='Depth of each scaffold of bin')
parser.add_argument("-o", '--output', required=True, help='Output file')

args = parser.parse_args()
## 
bin_file=open(args.bin)
contig_depth =open(args.depth)
outfile = open(args.output, 'w')

if not (os.path.exists(args.bin) or os.path.exists(args.depth)) :
    print "%s or %s was not exists, please check it."
    sys.exit()

#congtigs of bin dic 
bin_dic={}
for line in bin_file : 
    line_splt=line.strip().split("\t")
    if bin_dic.has_key(line_splt[1]):
        bin_dic[line_splt[1]].append(line_splt[0])
    else:
        bin_dic[line_splt[1]]=[line_splt[0]]

## load contig depth file
depth_dic={}
for contig in contig_depth : 
    contig_splt=contig.strip().split("\t")
    if contig_splt[0] == "contigName":
        continue
    depth_dic[contig_splt[0]]=contig_splt[1:3]


#calculate the weighted average depth
for k in bin_dic:
    numerator_total=0
    total_depth = 0
    total_length=0
    contig_names=bin_dic[k.strip()]
    for name in contig_names:
        try:
            contig_info=depth_dic[name]
        except KeyError:
            continue
        length=float(contig_info[0])
        depth=float(contig_info[1])
        numerator=length*depth
        numerator_total+=numerator
        total_length+=length
        total_depth+=depth

    try:
        average=numerator_total/total_length
    except ZeroDivisionError:
        continue
    out = [k, str(average), str(total_depth)]
    outfile.write("\t".join(out) + "\n")
