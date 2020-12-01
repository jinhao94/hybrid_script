#!/usr/bin/env python
'''
jinh, 2020
Merge files for Bin_cluster.py

Required files
1. Assemblied fasta file
2. depth file of assemblies
3. NT annotation info
4. UT annotation info
'''
import os, sys

if len(sys.argv) != 6:
    print "Usage: Bin_cluster_merge_files.py seq_fasta seq_depth nt_tax ut_tax out_prefix "
    sys.exit()

from Bio import SeqIO
import pandas as pd
from collections import defaultdict
from collections import Counter
import re

###### Functions
# def judge_nt_tax(id, c_tax):
#     if ';p_' in c_tax :
#         mc = re.search(r";p_(.*?);", c_tax)
#         p = mc.group(1)
#     else :
#         p = "None"
#     if ';c_' in c_tax :
#         mc = re.search(r';c_(.*?);', c_tax)
#         c = mc.group(1)
#     else :
#         c = "None"
#     if ';o_' in c_tax :
#         mc = re.search(r';o_(.*?);', c_tax)
#         o = mc.group(1)
#     else :
#         o = "None"
#     if ';f_' in c_tax :
#         mc = re.search(r';f_(.*?);', c_tax)
#         f = mc.group(1)
#     else :
#         f = "None"
#     if ';g_' in c_tax and id >= 80:
#         mc = re.search(r';g_(.*?);', c_tax)
#         g = mc.group(1)
#     else :
#         g = "None"
#     if ';s_' in c_tax and id >= 90:
#         mc = re.search(r';s_(.*)', c_tax)
#         s = mc.group(1)
#     else :
#         s = "None"
#     return [p, c, o, f, g, s]

def judge_nt_tax(id, total_len, aligned_len, c_tax):
    cov = float(aligned_len)/float(total_len)
    id = float(id)
    if cov < 0.5:
        return ["None", "None", "None", "None", "None", "None"]
    d,p,c,o,f,g,s,st = c_tax.split(";")
    tax_split = c_tax.split(";")
    if id >= 95 and "unclassified" not in tax_split:
        return [p, c, o, f, g, s]
    elif id >= 90 and g != "unclassified" :
        return [p, c, o, f, g, "None"]
    elif id >= 80 and f != "unclassified" :
        return [p, c, o, f, "None", "None"]
    elif id >=70 and o != "unclassified" :
        return [p, c, o, "None", "None", "None"]
    elif id >=60 and p != "unclassified" :
        return [p, "None", "None", "None", "None", "None"]
    else:
        return ["None", "None", "None", "None", "None", "None"]

'''
diamond file 
0 ~ 17
p pn pi c cn ci o on oi f fn fi g gn gi s sn si
'''

def judge_dia_tax(total_hits, c_tax):
    c_list = []
    if c_tax[0] != 'nan' and float(c_tax[1])/float(total_hits) > 0.5 :  # phylum
        p = c_tax[0]
    else: 
        p = "None"
    if c_tax[3] != 'nan' and float(c_tax[4])/float(total_hits) > 0.5 :  # class
        c = c_tax[3]
    else: 
        c = "None"
    if c_tax[6] != 'nan' and float(c_tax[7])/float(total_hits) > 0.5 : # order
        o = c_tax[6]
    else: 
        o = "None"
    if c_tax[9] != 'nan' and float(c_tax[10])/float(total_hits) > 0.5 and float(c_tax[11]) >= 60 : # family
        f = c_tax[9]
    else: 
        f = "None"
    if c_tax[12] != 'nan' and float(c_tax[13])/float(total_hits) > 0.5 and float(c_tax[14]) >= 70 : # genus
        g = c_tax[12]
    else: 
        g = "None"
    if c_tax[15] != 'nan' and float(c_tax[16])/float(total_hits) > 0.6 and float(c_tax[17]) >= 90 : # species
        s = c_tax[15]
    else: 
        s = "None"
    return [p, c, o, f, g, s]


fa = open(sys.argv[1], 'r')
dp = pd.read_csv(sys.argv[2], sep = "\t", header = None)
nt = pd.read_csv(sys.argv[3], sep = "\t", header = None)
ut = pd.read_csv(sys.argv[4], sep = "\t", header = None)

dp = dp.iloc[:, [0,2] ]
nt = nt.iloc[:, [0,1,4,5,6] ]

# set seq dict
seq = defaultdict(defaultdict)

for rec in SeqIO.parse(fa, 'fasta'):
    GC = (rec.seq.count("C") + rec.seq.count("G")) / float(len(rec.seq))
    seq[rec.description]['Length'] = len(rec.seq)
    seq[rec.description]['GC'] = GC

seq_feature = pd.DataFrame(seq).T
#seq_feature['seq_id'] = seq_feature.columns.values

seq_feature = pd.merge(seq_feature, dp, left_index =  True, right_on = 0, how="left") 
mg_out = pd.merge(seq_feature, nt, on = 0, how="left") 
mg_out = pd.merge(mg_out, ut, on = 0, how="left") 

# mg_out.to_csv("test.out", sep = "\t", index = True)
## for final tax
'''
col4 = nt id
col5 = nt tax
col6 = all proteins
col8 = get hit proteins
col9 = phylum
diamond file 
0 ~ 17
p pn pi c cn ci o on oi f fn fi g gn gi s sn si
'''

total_dict = defaultdict(defaultdict)
for index, row in mg_out.iterrows() :
    row_list = row.to_list()
    # print row_list
    seq_name = row_list[2]
    if str(row.iloc[6]) != 'nan' :
        total_dict[seq_name]['tax_list_nt'] = judge_nt_tax(row.iloc[6], row.iloc[1], row.iloc[5], row_list[7]) #id, total_len, aligned_len, c_tax
    else :
        total_dict[seq_name]['tax_list_nt'] = ["None"] * 6
    if str(row.iloc[10]) != 'nan' :
        total_dict[seq_name]['tax_list_ut'] = judge_dia_tax(row.iloc[10], row.to_list()[14:]) # total_hits, c_tax; 12: from phylum to species level
    else :
        total_dict[seq_name]['tax_list_ut'] = ["None"] * 6

seq_tax = pd.DataFrame(total_dict).T
# seq_tax.to_csv("test.tax", sep = "\t")

final_out = {}
for index, row in seq_tax.iterrows() :
    # print index
    if 'None' in row['tax_list_nt'] :
        if 'None' not in row['tax_list_ut'] :
            final_out[index] = row['tax_list_ut']
        else :
            nt = Counter(row['tax_list_nt'])
            ut = Counter(row['tax_list_ut'])
            if nt['None'] > ut['None'] :
                final_out[index] = row['tax_list_ut']
            else :
                final_out[index] = row['tax_list_nt']
    else :
         final_out[index] = row['tax_list_nt']


seq_feature.index = seq_feature[0]
seq_feature = seq_feature.drop(columns=[0])

outfile = pd.DataFrame(final_out).T
outfile = pd.merge(seq_feature, outfile, left_index=True, right_index=True, how = 'left')

outfile.to_csv(sys.argv[5], sep = "\t", na_rep = "None", header=False)
# seq_feature.to_csv("seq_feature.tax", sep = "\t")









