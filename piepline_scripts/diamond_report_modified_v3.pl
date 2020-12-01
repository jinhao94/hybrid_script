#!/usr/bin/env perl
#jinh 20191111
#update for getting clear notes for each contig/scaffolds or bins
################################################################
#可能存在一些bug, 有待反馈
#输出文件名问题, 待后续搭建流程时统一修改
#Future version: 提供最优注释信息
################################################################
use strict;
# use Bio::SeqIO;
use File::Basename;

sub help{
    print "Incorrect number of arguments\n Usage:  perl $0 bin.tax bin.prot.names bin.scaf2bin bin_name(or prefix) \n";
}

my $lineage_path="/gfsdata/gridengine/database/uniprot/uniprot_trembl.name_taxid_lineage";
unless (@ARGV==4) {
    help;
    exit;
}

my $tsv = shift;
my $Protein_path = shift;
my $Scaf2bin_path = shift;
my $prefix = shift;

# my $binid=$prefix;
$prefix=$prefix."_";

my %con;
my %bin;

# my $binid = basename $tsv;
my $obin=$prefix."_bin";
# print "Results will be written to $ocon and $obin\n";

# get gene number...
# my $in = Bio::SeqIO->new(-file => $faa, -format => 'fasta');
# while(my $seq = $in->next_seq) {
#     my $sid = $seq->primary_id;
#     my(@cdata) = split(/_/, $sid);
#     my $num = pop @cdata;
#     my $conid = join("_", @cdata);
#     # print "$conid\n";
#     $bin{$binid}{fap}++;
#     $con{$conid}{fap}++;
# }

my %scaf;

open(Scaf2bin, $Scaf2bin_path) or die "cannot open scaf2bin file... die";
while(<Scaf2bin>) {
    chomp();
    my @ss = split /\t/;
    $scaf{@ss[0]} = @ss[1];
}


open(Prot, $Protein_path) or die "cannot open protein list file... die";
while(<Prot>) {
    chomp();
    my @fs = split /\t/;
    my $conid = @fs[1];
    my $binid = $scaf{$conid};
    # print $conid;
    $bin{$binid}{fap}++;
    $con{$conid}{fap}++;
}


my %ids;
my $lastid = "none";

#parse lineage file
my %lineage;
my %h;

open(Lin, $lineage_path) or die "cannot open lineage files... die";
print "Parsing lineage database...\n";
while(<Lin>) {
    chomp();
    my @s = split /\t/;
    $h{@s[0]}=@s[1];
}

close Linn;
print "Getting results\n";

open(IN, $tsv) || die "cannot open tsv\n";
while(<IN>) {
    chomp();
    my ($qseqid,$sseqid,$stitle,$pident,$qlen,$slen,$length,$mismatch,$gapopen,$qstart,$qend,$sstart,$send,$evalue,$bitscore) = split(/\t+/);
    my $frac = $qlen / $slen;
    
    # print "\n\n\n$pident";
#get seq id  in diamond result file by delete _ conid means raw id (not be gene predicted id )
    my(@cdata) = split(/_/, $qseqid);
    my $id = pop @cdata;
    my $conid = join("_", @cdata);
    my $binid = $scaf{$conid};
    # print "DMD: $conid\n";

    if ($qseqid ne $lastid) {
        # we have a top hit

        my $org = undef;
        if ($stitle =~ m/OX=(\d+)\s/){
            $org = $1;
        } 
        # print "$org\n";
        my $tax = undef;
        if (exists $h{$org}){
            $tax = $h{$org};
        }
        my ($domain, $phylum, $class, $order, $family, $genus, $species) = split(/;/, $tax);

        $bin{$binid}{proteins}++;
        $con{$conid}{proteins}++;
        if ($frac > 0.7) {
            # print "$species\n";
            $bin{$binid}{fulllen}++;
            $con{$conid}{fulllen}++;
        #get count of different taxonomy
        # idendity
		# print "$domain\n"; 

        $bin{$binid}{phylum}{$phylum}{sump}+= $pident;
        $bin{$binid}{class}{$class}{sump}+= $pident;
        $bin{$binid}{order}{$order}{sump}+= $pident;
        $bin{$binid}{family}{$family}{sump}+= $pident;
        $bin{$binid}{genus}{$genus}{sump}+= $pident;
        $bin{$binid}{species}{$species}{sump}+= $pident;
        
        $bin{$binid}{phylum}{$phylum}{sumpn}++;
        $bin{$binid}{class}{$class}{sumpn}++;
        $bin{$binid}{order}{$order}{sumpn}++;
        $bin{$binid}{family}{$family}{sumpn}++;
        $bin{$binid}{genus}{$genus}{sumpn}++;
        $bin{$binid}{species}{$species}{sumpn}++;
        
        $con{$conid}{domain}{$domain}{sump}+= $pident;
        $con{$conid}{phylum}{$phylum}{sump}+= $pident;
        $con{$conid}{class}{$class}{sump}+= $pident;
        $con{$conid}{order}{$order}{sump}+= $pident;
        $con{$conid}{family}{$family}{sump}+= $pident;
        $con{$conid}{genus}{$genus}{sump}+= $pident;
        $con{$conid}{species}{$species}{sump}+= $pident;
        
		# print "$domain\t$con{$conid}{domain}{$domain}{sumpn}\n";
		$con{$conid}{domain}{$domain}{sumpn}++;
        $con{$conid}{phylum}{$phylum}{sumpn}++;
        $con{$conid}{class}{$class}{sumpn}++;
        $con{$conid}{order}{$order}{sumpn}++;
        $con{$conid}{family}{$family}{sumpn}++;
        $con{$conid}{genus}{$genus}{sumpn}++;
        $con{$conid}{species}{$species}{sumpn}++;

        # print "$conid\t$species\t$con{$conid}{species}{$species}{sumpn}\n"
    }
    # Top hist
    $lastid = $qseqid;
    }
}
close IN;

# unless (-d $outdir) {
    # mkdir $outdir;
# }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ for bin
open(BIN, ">$prefix".'bin');
while(my $bn = each(%bin)){
    # print $bn;
    print BIN $bn."\t".$bin{$bn}{fap}."\t".$bin{$bn}{proteins}."\t".$bin{$bn}{fulllen};
    
    my @sp = sort {$bin{$bn}{phylum}{$b}{sumpn} <=> $bin{$bn}{phylum}{$a}{sumpn} or $bin{$bn}{phylum}{$b}{sump} <=> $bin{$bn}{phylum}{$a}{sump}} keys %{$bin{$bn}{phylum}};
    my @sc = sort {$bin{$bn}{class}{$b}{sumpn} <=> $bin{$bn}{class}{$a}{sumpn} or $bin{$bn}{class}{$b}{sump} <=> $bin{$bn}{class}{$a}{sump}} keys %{$bin{$bn}{class}};
    my @so = sort {$bin{$bn}{order}{$b}{sumpn} <=> $bin{$bn}{order}{$a}{sumpn} or $bin{$bn}{order}{$b}{sump} <=> $bin{$bn}{order}{$a}{sump}} keys %{$bin{$bn}{order}};
    my @sf = sort {$bin{$bn}{family}{$b}{sumpn} <=> $bin{$bn}{family}{$a}{sumpn} or $bin{$bn}{family}{$b}{sump} <=> $bin{$bn}{family}{$a}{sump}} keys %{$bin{$bn}{family}};
    my @sg = sort {$bin{$bn}{genus}{$b}{sumpn} <=> $bin{$bn}{genus}{$a}{sumpn} or $bin{$bn}{genus}{$b}{sump} <=> $bin{$bn}{genus}{$a}{sump}} keys %{$bin{$bn}{genus}};
    my @ss = sort {$bin{$bn}{species}{$b}{sumpn} <=> $bin{$bn}{species}{$a}{sumpn} or $bin{$bn}{species}{$b}{sump} <=> $bin{$bn}{species}{$a}{sump}} keys %{$bin{$bn}{species}};


    #get mean idendity of each taxonomy, and print result of them.
    if ($bin{$bn}{phylum}{$sp[0]}{sumpn} > 0) {
        my $pmean = sprintf("%0.2f", $bin{$bn}{phylum}{$sp[0]}{sump} / $bin{$bn}{phylum}{$sp[0]}{sumpn});
        print BIN "\t", $sp[0], "\t", $bin{$bn}{phylum}{$sp[0]}{sumpn}, "\t", $pmean;
    } else {
        print BIN "\t", $sp[0], "\t", $bin{$bn}{phylum}{$sp[0]}{sumpn}, "\t", "0";
    }
    
    if ($bin{$bn}{class}{$sc[0]}{sumpn} > 0) {
        my $cmean = sprintf("%0.2f", $bin{$bn}{class}{$sc[0]}{sump} / $bin{$bn}{class}{$sc[0]}{sumpn});
        print BIN "\t", $sc[0], "\t", $bin{$bn}{class}{$sc[0]}{sumpn}, "\t", $cmean;
    } else {
        print BIN "\t", $sc[0], "\t", $bin{$bn}{class}{$sc[0]}{sumpn}, "\t", "0";
    }
    if ($bin{$bn}{order}{$so[0]}{sumpn} > 0) {
        my $omean = sprintf("%0.2f", $bin{$bn}{order}{$so[0]}{sump} / $bin{$bn}{order}{$so[0]}{sumpn});
        print BIN "\t", $so[0], "\t", $bin{$bn}{order}{$so[0]}{sumpn}, "\t", $omean;
    } else {
        print BIN "\t", $so[0], "\t", $bin{$bn}{order}{$so[0]}{sumpn}, "\t", "0";
    }
    
    if ($bin{$bn}{family}{$sf[0]}{sumpn} > 0) {
        my $fmean = sprintf("%0.2f", $bin{$bn}{family}{$sf[0]}{sump} / $bin{$bn}{family}{$sf[0]}{sumpn});
        print BIN "\t", $sf[0], "\t",$bin{$bn}{family}{$sf[0]}{sumpn}, "\t", $fmean;
    } else {
        print BIN "\t", $sf[0], "\t",$bin{$bn}{family}{$sf[0]}{sumpn}, "\t", "0";
    }
    
    if ($bin{$bn}{genus}{$sg[0]}{sumpn} > 0) {
        my $gmean = sprintf("%0.2f", $bin{$bn}{genus}{$sg[0]}{sump} / $bin{$bn}{genus}{$sg[0]}{sumpn});
        print BIN "\t", $sg[0], "\t", $bin{$bn}{genus}{$sg[0]}{sumpn}, "\t", $gmean;
    } else {
        print BIN "\t", $sg[0], "\t", $bin{$bn}{order}{$sg[0]}{sumpn}, "\t", "0";
    }
    
    if ($bin{$bn}{species}{$ss[0]}{sumpn} > 0) {
        my $smean = sprintf("%0.2f", $bin{$bn}{species}{$ss[0]}{sump} / $bin{$bn}{species}{$ss[0]}{sumpn});
        print BIN "\t", $ss[0], "\t", $bin{$bn}{species}{$ss[0]}{sumpn}, "\t", $smean;
    } else {
        print BIN "\t", $ss[0], "\t", $bin{$bn}{species}{$ss[0]}{sumpn}, "\t", "0";
    }
    print BIN "\n";
}
close BIN;

#exit;

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ for contigs/scaffolds
open(CON, ">$prefix".'con');
while(my($cn,$hr) = each %con) {
    # print "$cn\tdada\n"; 
    my $fl=$con{$cn}{fulllen};
    if ( $fl > 0 ) {
        print CON $cn."\t".$con{$cn}{fap}."\t".$con{$cn}{proteins}."\t".$con{$cn}{fulllen};
        my @sd = sort {$con{$cn}{domain}{$b}{sumpn} <=> $con{$cn}{domain}{$a}{sumpn} or $con{$cn}{domain}{$b}{sump} <=> $con{$cn}{domain}{$a}{sump}} keys %{$con{$cn}{domain}};
        my @sp = sort {$con{$cn}{phylum}{$b}{sumpn} <=> $con{$cn}{phylum}{$a}{sumpn} or $con{$cn}{phylum}{$b}{sump} <=> $con{$cn}{phylum}{$a}{sump}} keys %{$con{$cn}{phylum}};
        my @sc = sort {$con{$cn}{class}{$b}{sumpn} <=> $con{$cn}{class}{$a}{sumpn} or $con{$cn}{class}{$b}{sump} <=> $con{$cn}{class}{$a}{sump}} keys %{$con{$cn}{class}};
        my @so = sort {$con{$cn}{order}{$b}{sumpn} <=> $con{$cn}{order}{$a}{sumpn} or $con{$cn}{order}{$b}{sump} <=> $con{$cn}{order}{$a}{sump}} keys %{$con{$cn}{order}};
        my @sf = sort {$con{$cn}{family}{$b}{sumpn} <=> $con{$cn}{family}{$a}{sumpn} or $con{$cn}{family}{$b}{sump} <=> $con{$cn}{family}{$a}{sump}} keys %{$con{$cn}{family}};
        my @sg = sort {$con{$cn}{genus}{$b}{sumpn} <=> $con{$cn}{genus}{$a}{sumpn} or $con{$cn}{genus}{$b}{sump} <=> $con{$cn}{genus}{$a}{sump}} keys %{$con{$cn}{genus} };
        my @ss = sort {$con{$cn}{species}{$b}{sumpn} <=> $con{$cn}{species}{$a}{sumpn} or $con{$cn}{species}{$b}{sump} <=> $con{$cn}{species}{$a}{sump}} keys %{$con{$cn}{species}};
        
        # print "$cn\t$sp[0]\t$con{$cn}{phylum}{$sp[0]}{sump}\n";
        
        #get mean idendity of each taxonomy, and print result of them.
		if ($con{$cn}{domain}{$sd[0]}{sumpn} > 0) {
            my $pmean = sprintf("%0.2f", $con{$cn}{domain}{$sd[0]}{sump} / $con{$cn}{domain}{$sd[0]}{sumpn});
            print CON "\t", $sd[0], "\t", $con{$cn}{domain}{$sd[0]}{sumpn}, "\t", $pmean;
        } else {
            print CON "\t", $sd[0], "\t", $con{$cn}{domain}{$sd[0]}{sumpn}, "\t", "0";
        }
		
        if ($con{$cn}{phylum}{$sp[0]}{sumpn} > 0) {
            my $pmean = sprintf("%0.2f", $con{$cn}{phylum}{$sp[0]}{sump} / $con{$cn}{phylum}{$sp[0]}{sumpn});
            print CON "\t", $sp[0], "\t", $con{$cn}{phylum}{$sp[0]}{sumpn}, "\t", $pmean;
        } else {
            print CON "\t", $sp[0], "\t", $con{$cn}{phylum}{$sp[0]}{sumpn}, "\t", "0";
        }
        
        if ($con{$cn}{class}{$sc[0]}{sumpn} > 0) {
            my $cmean = sprintf("%0.2f", $con{$cn}{class}{$sc[0]}{sump} / $con{$cn}{class}{$sc[0]}{sumpn});
            print CON "\t", $sc[0], "\t", $con{$cn}{class}{$sc[0]}{sumpn}, "\t", $cmean;
        } else {
            print CON "\t", $sc[0], "\t", $con{$cn}{class}{$sc[0]}{sumpn}, "\t", "0";
        }
        if ($con{$cn}{order}{$so[0]}{sumpn} > 0) {
            my $omean = sprintf("%0.2f", $con{$cn}{order}{$so[0]}{sump} / $con{$cn}{order}{$so[0]}{sumpn});
            print CON "\t", $so[0], "\t", $con{$cn}{order}{$so[0]}{sumpn}, "\t", $omean;
        } else {
            print CON "\t", $so[0], "\t", $con{$cn}{order}{$so[0]}{sumpn}, "\t", "0";
        }
        
        if ($con{$cn}{family}{$sf[0]}{sumpn} > 0) {
            my $fmean = sprintf("%0.2f", $con{$cn}{family}{$sf[0]}{sump} / $con{$cn}{family}{$sf[0]}{sumpn});
            print CON "\t", $sf[0], "\t", $con{$cn}{family}{$sf[0]}{sumpn}, "\t", $fmean;
        } else {
            print CON "\t", $sf[0], "\t", $con{$cn}{family}{$sf[0]}{sumpn}, "\t", "0";
        }
        
        if ($con{$cn}{genus}{$sg[0]}{sumpn} > 0) {
            my $gmean = sprintf("%0.2f", $con{$cn}{genus}{$sg[0]}{sump} / $con{$cn}{genus}{$sg[0]}{sumpn});
            print CON "\t", $sg[0], "\t", $con{$cn}{genus}{$sg[0]}{sumpn}, "\t", $gmean;
        } else {
            print CON "\t", $sg[0], "\t", $con{$cn}{order}{$sg[0]}{sumpn}, "\t", "0";
        }
        
        if ($con{$cn}{species}{$ss[0]}{sumpn} > 0) {
            my $smean = sprintf("%0.2f", $con{$cn}{species}{$ss[0]}{sump} / $con{$cn}{species}{$ss[0]}{sumpn});
            print CON "\t", $ss[0], "\t", $con{$cn}{species}{$ss[0]}{sumpn}, "\t", $smean;
        } else {
            print CON "\t", $ss[0], "\t", $con{$cn}{species}{$ss[0]}{sumpn}, "\t", "0";
        }
        print CON "\n";
    }
}

close CON;
