# hybrid_script
Script for hybrid
# Hybrid maintext scripts
Hybrid, extra-deep metagenomic sequencing enables genomic and functional characterization of low-abundance and uncultured species in the human gut microbiome

# Hybrid assembly and binning pipeline
Getting Started 

## Dependencies
You will need the following python packages installed.

[Biopython](https://biopython.org/)

[Pandas](https://pandas.pydata.org/)

**Software**

[Diamond](https://github.com/bbuchfink/diamond)

[BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastNews)

[Flye 2.7](https://github.com/fenderglass/Flye)

[SPAdes 3.14](https://github.com/ablab/spades)

[NextPolish](https://github.com/Nextomics/NextPolish)

[Quickmerge](https://github.com/mahulchak/quickmerge)

[Minimap](https://github.com/lh3/minimap)

[Bowtie2](https://github.com/BenLangmead/bowtie)

[berokka](https://github.com/tseemann/berokka)

[csvtk](https://github.com/shenwei356/csvtk)

## Preprocessing
## Step1 Assembly
```
# Run flye
flye --meta --pacbio-raw Pacbio.fasta -m 300m -t 32 -o 01_Pacbio_assembly
# Polish Pacbio_assemblies (nextpolish.sh is in Assembly_script)
nextpolish.sh -r1 shrot_read_r1.fq -r2 shrot_read_r2.fq -lr Pacbio.fasta -a 01_Pacbio_assembly/assembly.fasta -o 02_Pacbio_assembly_polished
nextPolish 02_Pacbio_assembly_polished/run.cfg
# Run SPAdes
SPAdes --pacbio Pacbio.fasta -m 300m -t 32 -o Pacbio_assembly -1 shrot_read_r1.fq -2 shrot_read_r2.fq -t 32 -o 03_hybrid_out
# Merge assembles using quickmerge (merge_assembles.sh is in Assembly_script)
merge_assembles.sh -p 03_hybrid_out/scaffolds.fasta -s 02_Pacbio_assembly_polished/genome.nextpolish.fasta -o 04_merged_assembly -t 8
```
## Step2 Error correction
The assembly error correction step aimed at localizing inconsistent coverage region and cut off the sequence at breakpoint.
1. Clear overhang in assemblies
```
berokka --outdir 02_trim_overhang/02.trimmed.fa 04_merged_assembly/merged.fasta
```
2. Caculate the region depth of short and long read.
```
mkdir 05_region_depth
minimap2 -t 32 -aLQx map-pb --secondary=no --sam-hit-only 02_trim_overhang/02.trimmed.fa Pacbio.fasta | samtoolsview -bS -@ 10 - | samtools sort - -@ 5 -m 10g -o 05_region_depth/long.sort.bam
bowtie-build 02_trim_overhang/02.trimmed.fa -p 02_trim_overhang/02.trimmed.fa
bowtie2 --no-unal --no-discordant --end-to-end --very-sensitive -x 02_trim_overhang/02.trimmed.fa -1 shrot_read_r1.fq -2 shrot_read_r1.fq -p 60 | samtools view -bS -@ 16 - | samtools sort -@ 16 -m3g - -o  05_region_depth/short.sort.bam

mosdepth -t 32 -b 1000 -n -F 256 05_region_depth/long.depth 05_region_depth/long.sort.bam
mosdepth -t 32 -b 1000 -i 2 -F 1796 05_region_depth/short.depth 05_region_depth/short.sort.bam

gunzip 05_region_depth/*depth.regions.bed.gz
for i in `ls -d 05_region_depth/*depth.regions.bed`; less $i | perl -lane 'print "@F[0]_@F[1]\t$_"' > ${i}.f ; done 

csvtk join -t -T -H -k 05_region_depth/short.depth.regions.bed.f 05_region_depth/long.depth.regions.bed.f > depth.merged

# Gene prediction
prodigal -i 02_trim_overhang/02.trimmed.fa -a 02_trim_overhang/02.trimmed.faa -f gff -o 02_trim_overhang/02.trimmed.gff
```

```
python Assembly_script/correct_scaf.py 
usage: correct_scaf.py [-h] -t TABLE -g GFF -a FASTA -o OUTDIR [-f] [-l LOWER]
                       [-u UPPER] [-i INTERVAL] [-m MINLEN]
correct_scaf.py: error: argument -t/--table is required

python Assembly_script/correct_scaf.py  -t depth.merged -g 02_trim_overhang/02.trimmed.gff -a 02_trim_overhang/02.trimmed.fasta -o merged.corrected/corrected.fa
``` 

3. filt the sequence length lower than 2Kbp
```

18_select_target_length_seq.py merged.corrected/corrected.fa merged.corrected/corrected_minlen2K.fa

```
## Step3 Binning
1. Run metabat2 
```
bwa index merged.corrected/corrected_minlen2K.fa 
bwa mem -t 32 merged.corrected/corrected_minlen2K.fa shrot_read_r1.fq shrot_read_r2.fq | samtools view -bS - | samtools sort -@ 16 - -o corrected_minlen2K.sort.bam
jgi_summarize_bam_contig_depths --outputDepth corrected_minlen2K.sort.depth corrected_minlen2K.sort.bam
mkdir 03_binning/metabat
metabat2 -i merged.corrected/corrected_minlen2K.fa -a corrected_minlen2K.sort.depth -o 03_binning/metabat/metabat --minContig 2000 --unbinned

```
2. Run Cluster
This step required a per-bulit UniProt TrEMBL and NCBI nt database 
```
mkdir 03_binning/nt_ut_tax
# predicted the coding genes of contigs
prodigal -i merged.corrected/corrected_minlen2K.fa -f gff -a merged.corrected/corrected_minlen2K.faa -q -o merged.corrected/corrected_minlen2K.gff 
# Run diamond
diamond blastp --threads 100 --max-target-seqs 10 --db  uniprot_trembl --query merged.corrected/corrected_minlen2K.faa --outfmt 6 qseqid sseqid stitle pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore --out 03_binning/nt_ut_tax/corrected_minlen2K.tsv
diamond_report_modified_v2.pl 03_binning/nt_ut_tax/corrected_minlen2K.tsv merged.corrected/corrected_minlen2K.faa 03_binning/nt_ut_tax/corrected_minlen2K.ut_tax

# Run blastn
Run_blast_piep.py -i merged.corrected/corrected_minlen2K.fa -o 03_binning/nt_ut_tax/corrected_minlen2K.nt_tax -t 32 -p corrected_minlen2K 

# extract single-copy core genes for each contigs
mkdir 03_binning/nt_ut_tax/SGC
hmmscan -o 03_binning/nt_ut_tax/SGC/corrected_minlen2K.hmm.out --tblout 03_binning/nt_ut_tax/SGC/corrected_minlen2K.blast.out --cpu 6 /nvmessdnode3/opt/database/123_SCG_hmm_db/bacteria_139_CSCG.hmm merged.corrected/corrected_minlen2K.faa

# calculate the 5mer frequency 
create_kmer_freq_modified.py merged.corrected/corrected_minlen2K.fa 5 > 03_binning/nt_ut_tax/corrected_minlen2K.5mer
## step1 Input files preparation 
Bin_cluster_merge_files.py merged.corrected/corrected_minlen2K.fa corrected_minlen2K.sort.depth 03_binning/nt_ut_tax/corrected_minlen2K.nt_tax/corrected_minlen2K.anno.f 03_binning/nt_ut_tax/corrected_minlen2K.ut_tax_con 03_binning/cluster_file.in

## step2 Run contig cluster
Bin_cluster.py 03_binning/cluster_file.in 03_binning/nt_ut_tax/corrected_minlen2K.5mer 03_binning/nt_ut_tax/SGC/corrected_minlen2K.blast.out 03_binning/cluster_out

## step3 Run checkM
checkm lineage_wf -x fa -t 32 03_binning/metabat 03_binning/metabat_checkm
checkm lineage_wf -x fa -t 32 03_binning/cluster_out 03_binning/cluster_out_checkm 

summarize_checkm.py 03_binning/metabat_checkm > 03_binning/metabat_checkm.s 
summarize_checkm.py 03_binning/cluster_out_checkm > 03_binning/cluster_out_checkm.s 

## step4 choose best bins 
Bin_choose_best.py -a merged.corrected/corrected_minlen2K.fa -c 03_binning/cluster_out -m 03_binning/metabat/metabat -cc 03_binning/cluster_out_checkm.s -mc 03_binning/metabat_checkm.s -o 03_binning/integrated_bins
checkm lineage_wf -x fa -t 32 03_binning/integrated_bins 03_binning/integrated_bins_checkm 
summarize_checkm.py 03_binning/integrated_bins_checkm > 03_binning/integrated_bins_checkm.s
sh extract_faa_name.sh merged.corrected/corrected_minlen2K.faa '>' merged.corrected/corrected_minlen2K.prot.names

## step5 Merge bins
Bin_combined_subbins.py -r 03_binning/integrated_bins -t 03_binning/nt_ut_tax/corrected_minlen2K.ut_tax_con -d corrected_minlen2K.sort.depth -pt merged.corrected/corrected_minlen2K.prot.names -c 03_binning/integrated_binsintegrated_bins_checkm.s -o 03_binning/merge_bin -p merge_bin

checkm lineage_wf -x fa -t 32 03_binning/merge_bin 03_binning/merge_bin_checkm
summarize_checkm.py 03_binning/merge_bin_checkm > 03_binning/merge_bin_checkm.s
## step6 choose best bins 
bestbin.sh 03_binning/merge_bin_checkm.s 03_binning/merge_bin/merge_bin.info 03_binning/merge_bin/merge_bin.mag_file 03_binning/integrated_bins 03_binning/merge_bin/merge_bin_combine_bins 03_binning/final_bin
