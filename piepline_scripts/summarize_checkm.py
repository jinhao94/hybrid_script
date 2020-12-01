#!/usr/bin/env python
import sys
import os
import logging
# This script summarizes the statistics of each bin by parsing 
# the checkm_folder/storage/bin_stats_ext.tsv file of the CheckM output
if len(sys.argv) != 2 :
    print "Usage: summarize_checkm.py checkm_out_folder" 
    sys.exit()

bin_file = sys.argv[1] + "/storage/bin_stats_ext.tsv"
if not os.path.exists(bin_file):
    logging.error("bin_stats_ext.tsv not exists, please checkm your input checkm path or checkm result")
    sys.exit()

for line in open(bin_file):
    dic=eval(line.strip().split("\t")[1])
    #if dic["Completeness"]<20 or dic["Contamination"]>10: continue
    if "__" in dic["marker lineage"]: dic["marker lineage"]=dic["marker lineage"].split("__")[1]
    #if "unbinned" in line.split("\t")[0]: name="Unbinned"
    name=line.split("\t")[0]
    print "\t".join([name, str(dic["Completeness"])[:5],\
    str(dic["Contamination"])[:5], str(dic["GC"])[:5],\
    dic["marker lineage"], str(dic["N50 (contigs)"]),\
    str(dic["Genome size"])])
