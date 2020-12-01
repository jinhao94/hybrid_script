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


parser = argparse.ArgumentParser(description='Merge the incompleteness bins.')
parser.add_argument("-r", '--rawbins', required=True, help='rawbins fasta')
parser.add_argument("-t", '--tax', required=True, help='taxonomy of contig')
parser.add_argument("-d", '--depth', required=True, help='depth of contig')
parser.add_argument("-pt", '--protein', required=True, help='protein names file of contig')
parser.add_argument("-c", '--checkm', required=True, help='checkm out')
parser.add_argument("-o", '--outdir', required=True, help='output file')
parser.add_argument("-p", '--prefix', required=True, help='output prefix(sample name)')
parser.add_argument("-f", '--force', action = 'store_true', help='overwriting the existing files')

args = parser.parse_args()

rawbins = args.rawbins
tax = args.tax
depth = args.depth
protein = args.protein
outdir = args.outdir
prefix = args.prefix
checkm = args.checkm

if not os.path.exists(outdir):
    os.makedirs(outdir)
else:
    if args.force == True:
        shutil.rmtree(outdir)
        os.makedirs(outdir)
    else:
        print "Outdir already exist, please rename or use --force"
        sys.exit()

## Required files processing
scaf2bin = outdir + "/" + prefix + ".scaf2bin" ## s1 scaf2bin
bin_depth = outdir + "/" + prefix + ".bin_depth" ## s2 depth 
tax_prefix = outdir + "/" + prefix  ## s3 taxmonmy 
tax_bin = outdir + "/" + prefix  + "_bin" ## s3 bin tax
bin_info = outdir + "/" + prefix + ".info"  ## s4 bin_info
bin_outprefix = outdir + "/" + prefix
# bin_cls = outdir + "/" + prefix + ".cls"  ## s5 bin_cls
# bin_mag_list = outdir + "/" + prefix + ".mag.list"  ## s6 bin_cls
# bin_cls_list = outdir + "/" + prefix + ".cls.list" ## s7 bin_cls.list
# bin_mag_file = outdir + "/" + prefix + ".mag.file" ## s8 bin_mag_file
# bin_outdir = outdir + "/" + prefix + ".merged_bins"


## Commond
scaf2bin_cmd  = "Fasta_to_Scaffolds2Bin.sh -i %s > %s" %(rawbins, scaf2bin)
bin_depth_cmd = "Bin_depth.py -b %s -d %s -o %s "  %(scaf2bin, depth, bin_depth)
tax_cmd = "diamond_report_modified_v3.pl %s %s %s %s" %(tax, protein, scaf2bin, tax_prefix)
bin_info_cmd = 'less %s | cut -f1,2,3,4,6,7 | csvtk join -t -T -H - %s | csvtk join -t -T -H - %s | cut -f1-7,11,27-29 | sort -k9,9 -k7,7nr -k2,2nr -k3,3n > %s' %(checkm, bin_depth, tax_bin, bin_info)

## Run merge bin pipeline
merge_cmd = "mergebin.sh %s %s %s" %(bin_info, rawbins, bin_outprefix)
# s5_cmd = 'bin.clust2.pl %s %s' %(bin_info, bin_cls)
# s6_cmd = ''' grep -P '^CCC' Sample_LCB10.int.cls | perl -ne  'chomp;@s=split /\t+/; shift @s;$s[1]=100 if $s[1]>100; print "".(join "\t",@s)."\n";' > %s ''' %(bin_cls, bin_mag_list)
# s7_cmd = """ grep -vP "^FIL" Sample_LCB10.int.cls | perl -ne 'chomp;@s=split /\s+/;print "$s[1]\t" if $s[0] eq "CON";print "$s[1]\n" if $s[0] eq "CCC";' | perl -ne'chomp;@s=split /\s+/;$a=pop @s;print "$a\t".($#s+1)."\t".(join ";",@s)."\n";' > %s """ %(bin_cls, bin_cls_list)

# s8_cmd = "diff_py %s 1 %s 1 %s" %(bin_mag_list, bin_cls_list, bin_mag_file)

# s9_cmd = """ mkdir %s; perl -e 'open I, @ARGV[0]; while(<I>){chomp; @s=split /\s+/; $o=@ARGV[2]."/".@s[0]; @b=split /;/, @s[2]; $o="cat @ARGV[1]/".(join ".fa @ARGV[1]/", @b).".fa > @ARGV[2]/@s[0].fa"; print "$o\n"}' %s %s %s """ %(bin_outdir, bin_mag_file, rawbins, bin_outdir)

for i in [scaf2bin_cmd, bin_depth_cmd, tax_cmd, bin_info_cmd, merge_cmd]:
    subprocess.call(i, shell = True)




