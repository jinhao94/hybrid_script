#!/bin/sh
if [ $# -ne 6 ];then
    echo "Usage bestbin.sh merged_bin_checkm bin_info magfile raw_bins merged_bins outdir"
    exit
fi
merged_bin_checkm=$1
bin_info=$2
magfile=$3
raw_bins=$4
merged_bins=$5
outdir=$6
mkdir $outdir
perl -e 'open CCM, @ARGV[0]; while(<CCM>){chomp; @s=split /\s+/; $ccm_score = @s[1]-5*@s[2]; $ccm{@s[0]}{"score"}=$ccm_score;}; open RCM, @ARGV[1]; while(<RCM>){chomp; @s=split /\s+/; $rcm_score = @s[1]-5*@s[2]; $rcm{@s[0]}{"score"}=$rcm_score;}; open MAG, @ARGV[2]; while(<MAG>){chomp; @s=split /\t/;if(not $_=~/;/){ print "@s[0]\n"; next} ; @s1 = split /;/, @s[2]; foreach $a (@s1){$cs=$ccm{@s[0]}{"score"}; $rs=$rcm{$a}{"score"}; if($rs >= $cs){$o=$a; break}else{$o=@s[0]} }; print "$o\n"}' $1 $2 $3 | grep 'cbin' | parallel -j 5 cp ${merged_bins}/{}.fa $outdir

perl -e 'open CCM, @ARGV[0]; while(<CCM>){chomp; @s=split /\s+/; $ccm_score = @s[1]-5*@s[2]; $ccm{@s[0]}{"score"}=$ccm_score;}; open RCM, @ARGV[1]; while(<RCM>){chomp; @s=split /\s+/; $rcm_score = @s[1]-5*@s[2]; $rcm{@s[0]}{"score"}=$rcm_score;}; open MAG, @ARGV[2]; while(<MAG>){chomp; @s=split /\t/;if(not $_=~/;/){ print "@s[0]\n"; next} ; @s1 = split /;/, @s[2]; foreach $a (@s1){$cs=$ccm{@s[0]}{"score"}; $rs=$rcm{$a}{"score"}; if($rs >= $cs){$o=$a; break}else{$o=@s[0]} }; print "$o\n"}' $1 $2 $3 | grep -v 'cbin' | parallel -j 5 cp ${raw_bins}/{}.fa $outdir 

