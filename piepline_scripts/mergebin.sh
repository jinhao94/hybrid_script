#!/bin/sh
bin_file=$1
raw_bin_fa=$2
outprefix=$3
 
sort -k9,9 -k7,7nr -k2,2nr -k3,3n $1 -o $1

bin.clust2.pl $1 ${3}.cls

grep -P '^CCC' ${3}.cls | perl -ne  'chomp;@s=split /\t+/; shift @s;$s[1]=100 if $s[1]>100;print "".(join "\t",@s)."\n";' > ${3}.mag.list

grep -vP "^FIL" ${3}.cls | perl -ne 'chomp;@s=split /\s+/;print "$s[1]\t" if $s[0] eq "CON";print "$s[1]\n" if $s[0] eq "CCC";' | perl -ne 'chomp;@s=split /\s+/;$a=pop @s;print "$a\t".($#s+1)."\t".(join ";",@s)."\n";' > ${3}.cls.list

diff_py ${3}.mag.list 1 ${3}.cls.list 1 ${3}.mag_file

mkdir ${3}_combine_bins; perl -e 'open I, @ARGV[0]; while(<I>){chomp; @s=split /\s+/; $o=@ARGV[2]."/".@s[0]; @b=split /;/, @s[2]; $o="cat @ARGV[1]/".(join ".fa @ARGV[1]/", @b).".fa > @ARGV[2]/@s[0].fa"; print "$o\n"}' ${3}.mag_file $2 ${3}_combine_bins | parallel {}
