#!/usr/bin/env python
## jinh
import os
import sys
from collections import defaultdict
from Bio import SeqIO
import argparse
import shutil
import subprocess

parser = argparse.ArgumentParser(description="A pipeline for aligning the sequence with NCBI-nt database. \n Use -h or -help to print detailed descriptions of command line arguments ")

if len(sys.argv)==1:
    print parser.print_help()

parser.add_argument("-i", '--fasta', required=True, help='Fasta file to be aligned against with nt database')
parser.add_argument("-o", '--outdir', required=True, help='Output directory')
# parser.add_argument("-ot", '--out_type', default="species", choices=["sub", "subspecies", "species", "sp"] , help='Output type. Get annotation at species or subspecies level')
parser.add_argument("-f", '--force' , default= False ,action = 'store_true', help='Overwriting the existing files')
parser.add_argument("-t", '--threads' , default= 16 , type = int, help='Number of threads (default 16)')
parser.add_argument("-p", '--prefix' , default= 'blast' , type = str, help='Prefix (default blast)')
parser.add_argument("-d", '--dryrun' , default= False ,action = 'store_true', help='Only print command')

args = parser.parse_args()

fasta = args.fasta
# out_type = args.out_type
threads = args.threads
froce = args.force
outdir = args.outdir
prefix = args.prefix
dryrun = args.dryrun

blast_out = outdir + "/" + prefix + ".bt"
blast_conn = blast_out + ".conn"

if not os.path.exists(outdir):
    os.makedirs(outdir)
else:
    if args.force == True:
        shutil.rmtree(outdir)
        os.makedirs(outdir)
    else:
        print "Outdir already exist, please rename or use --force"
        sys.exit()


# table_tmp = open(outdir + "/table.tmp", 'w')
# table_tmp_all = open(outdir + "/table.tmp_all", 'w')
# readme_file = "/nvmessdnode3/opt/.method/.script/NT_uniport/get_species.anno.readme"

## run blastn
blast_cmd = "blastn -query %s -db //gfsdata/gridengine/database/nt_v5/nt -out %s -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid' -evalue 1e-10 -word_size 24 -max_target_seqs 10 -num_threads %s" %(fasta, blast_out, threads)
## connect hsp
conn_cmd = 'conn_blast.py %s %s' %(blast_out, blast_conn)
## seq_len
seq_len = outdir + "/" + fasta + ".len"
blast_anno = outdir + "/" + fasta + ".anno"
seq_len_cmd = 'seq_len %s > %s' %(fasta, seq_len)
get_species_cmd = 'get_species.anno.py -b %s -l %s -o %s' %(blast_conn, seq_len, blast_anno)


if dryrun == True:
    print blast_cmd
    print conn_cmd
    print seq_len_cmd
    print get_species_cmd
else:
    subprocess.call(blast_cmd, shell=True)
    subprocess.call(conn_cmd, shell=True)
    subprocess.call(seq_len_cmd, shell=True)
    subprocess.call(get_species_cmd, shell=True)





