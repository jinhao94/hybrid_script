#!/usr/bin/env python
import os
import sys
from collections import defaultdict
from Bio import SeqIO
import argparse
import shutil
import subprocess
import cPickle
import logging

# gb_ass_file = "/nvmessdnode3/opt/database/taxdump/nucl_gb.accession2taxid"
# tax_species = "/nvmessdnode3/opt/database/taxdump/nucl_gb.accession2taxid.species"
tax_file_path = "/nvmessdnode3/opt/database/taxdump/tax_id.strain.lineage"

parser = argparse.ArgumentParser(description='Get annotation at species/subspecies level.')

parser.add_argument("-b", '--blast', required=True, help='blast connected table')
# parser.add_argument("-t", '--type', default= "species", choices=["sub", "subspecies", "species", "sp"] , help='get annotation at species or subspecies level (default = species)')
parser.add_argument("-l", '--length', required=True, help='sequence length file (A file containing the name and length of the sequence and separated by TAB)')
parser.add_argument("-o", '--outfile', required=True, help='output file')
parser.add_argument("-f", '--force' , default= False ,action = 'store_true', help='overwriting the existing files')


args = parser.parse_args()

blast = args.blast
# out_type = args.type
lensfile = args.length
outfile = args.outfile
force = args.force

## check files
if not os.path.exists(tax_file_path):
    logging.error('Core file not exists, please check %s.' %(tax_file_path)) 
    sys.exit()


# sort by score

seq_len = {}
for seq in open(lensfile):
    seq_sp = seq.strip().split("\t")
    seq_len[seq_sp[0]] = seq_sp[1]

## load or create tax dict
tax_pkl = tax_file_path + ".pkl"
if not os.path.exists(tax_pkl):
    tax_file = open(tax_file_path)
    tax_dic = {}
    for line in tax_file:
        line_sp = line.strip().split("\t")
        tax_dic[line_sp[0]] = line_sp[1:]
    cPickle.dump(tax_dic, open(tax_pkl, "wb"))
else:
    tax_dic = cPickle.load(open(tax_pkl))

# print "tax_dic done"

out_dic = defaultdict(dict)

for line in open(blast):
    line_sp = line.strip().split("\t")
    contig, gene, idendity, align_len= line_sp[0:4]
    tax_id = line_sp[9].split(";")[0]
    try:
        tax_species = tax_dic[tax_id]
    except KeyError:
        try:
            tax_id_sup = line_sp[9].split(";")[1]
        except IndexError:
            continue
        tax_species =tax_dic[tax_id_sup]
    tax=tax_species[1].split(";")
    specie_taxid = tax[6]
    strain_taxid = tax[7]
    if contig not in out_dic:
        out_dic[contig]["total"] = 1
        out_dic[contig]["species_counter"] = 1
        out_dic[contig]["strain_counter"] = 1
        out_dic[contig]["strain_total"] = 0
        out_dic[contig]["species_total"] = 0
        proportion = str(out_dic[contig]["strain_counter"])  + "/" + str(out_dic[contig]["total"]) + "|" + str(out_dic[contig]["species_counter"])  + "/" + str(out_dic[contig]["total"]) 
        out_dic[contig]["out"] = [contig, seq_len[contig], gene, proportion, align_len, idendity, ";".join(tax)]
        continue
    else:
        # print tax
        ctg_species_taxid = out_dic[contig]['out'][6].split(";")[6]
        ctg_strain_taxid = out_dic[contig]['out'][6].split(";")[7]
        if  strain_taxid != "unclassified" and strain_taxid == ctg_strain_taxid:
            out_dic[contig]["strain_counter"] += 1
        else:
            out_dic[contig]["strain_total"] +=1
        if specie_taxid != "unclassified" and specie_taxid == ctg_species_taxid:
            out_dic[contig]["species_counter"] += 1
        else:
            out_dic[contig]["species_total"] +=1
    out_dic[contig]["total"] += 1
    proportion = str(out_dic[contig]["strain_counter"])  + "/" + str(out_dic[contig]["strain_total"]) + "|" + str(out_dic[contig]["species_counter"])  + "/" + str(out_dic[contig]["species_total"]) + "|" + str(out_dic[contig]["total"])
    out_dic[contig]["out"][3] = proportion


output = open(outfile, 'w')
header = ['seq', 'length', 'gene_id', 'same_top_strain/total_strain|same_top_species/total_species|total_hits', 'align_length', 'idendity', 'lineage']
output.write("\t".join(header) + "\n")
for key in sorted(out_dic.keys()):
    out_join =  "\t".join(out_dic[key]['out'])
    output.write(out_join + "\n")


