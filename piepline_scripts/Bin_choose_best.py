#!/usr/bin/env python
import os
import sys
from Bio import SeqIO
import argparse
import shutil
import subprocess
import logging
from collections import defaultdict
import subprocess


parser = argparse.ArgumentParser(description='Choose best bin for contigs')
parser.add_argument("-a", '--assembly_fasta', required=True, help='assembled fasta')
parser.add_argument("-c", '--cluster', required=True, help='cluster bins ')
parser.add_argument("-m", '--metabat', required=True, help='metabat bins')
parser.add_argument("-cc", '--cluster_cm', required=True, help='checkm result of cluster_method')
parser.add_argument("-mc", '--metabat_cm', required=True, help='checkm result of metabat2 bins')
parser.add_argument("-o", '--outdir', required=True, help='output file')
parser.add_argument("-f", '--force', action = 'store_true', help='overwriting the existing files')

args = parser.parse_args()

assembly = args.assembly_fasta
cluster_cm = args.cluster_cm
metabat_cm = args.metabat_cm
cluster = args.cluster
metabat = args.metabat
outdir = args.outdir

## Check status of files
if not (os.path.exists(assembly) or os.path.exists(cluster) or os.path.exists(metabat) ):
    logging.error('Core file not exit, please check %s, %s.' %(assembly, binfile)) 
    sys.exit()

if not os.path.exists(outdir):
    os.makedirs(outdir)
else:
    if args.force == True:
        shutil.rmtree(outdir)
        os.makedirs(outdir)
    else:
        print "Outdir already exist, please rename or use --force"
        sys.exit()

## Merge files
# Output files
scaf2bin_mb = outdir + "/" + "metabat.scaf2bin"
scaf2bin_c = outdir + "/" + "cluster.scaf2bin"
binfile = outdir + "/" + "merged_file"

# Run commands
scaf2bin_mb_cmd = "Fasta_to_Scaffolds2Bin.sh -i %s > %s" %(metabat, scaf2bin_mb)
scaf2bin_c_cmd = "Fasta_to_Scaffolds2Bin.sh -i %s > %s" %(cluster, scaf2bin_c)
merge_file_cmd = 'csvtk join -t -T -H -k %s %s | csvtk join -t -T -H -k -f"2;1" - %s | csvtk join -t -T -H -k -f"3;1" - %s > %s' %(scaf2bin_c, scaf2bin_mb, cluster_cm, metabat_cm, binfile)

subprocess.call(scaf2bin_mb_cmd, shell=True)
subprocess.call(scaf2bin_c_cmd, shell=True)
subprocess.call(merge_file_cmd, shell=True)


binfile_read = open(binfile, 'r')
out_list  = []
for line in binfile_read:
    line_sp = line.strip().split("\t")
    Ccon, Ccom = line_sp[3:5]
    try :
        Mcon, Mcom = line_sp[9:11]
    except ValueError:
        Mcon, Mcom = [100,10000]
    Cs = float(Ccon) - 5*float(Ccom)
    Ms = float(Mcon) - 5*float(Mcom)
    if Cs >= Ms:
        Best = line_sp[1]
    else:
        Best = line_sp[2]
    out_list.append([Best] + line_sp)

# print out_list

out_file_name = outdir + "/bin_best.txt"
outfa_dir = outdir + "/best_fasta"
os.makedirs(outfa_dir)
bin_dic = defaultdict()
out_file = open(out_file_name, 'w')

for i in out_list:
    bin_dic[i[1]] = i[0]
    out = "\t".join(i)
    out_file.write(out + "\n")

# print bin_dic

for seq in SeqIO.parse(assembly, "fasta"):
    if seq.id in bin_dic:
        outfile = outfa_dir + "/" + bin_dic[seq.id] + '.fa'
        if not os.path.exists(outfile):
            outfile_wt = open(outfile, 'w')
        else:
            outfile_wt = open(outfile, 'aw')
        outfile_wt.write(">" + str(seq.id) + "\n" + str(seq.seq) + "\n")




