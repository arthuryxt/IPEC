#!/usr/bin/perl -w
use strict;
# From SMRTanalysis 2 onward, No raw sequences are provided. Only subreads and ccs reads are provided.
# There will be three movies: p0.1, p0.2 and p0.3
# There will be only s1 (as the CCD is recoding ALL the ZMW ALL the time)
# report file-1 : ccs-reads
# report file-2 : subreads with NO ccs-reads
die "Usage: $0  \" dir\"  \"serial\/lane number\"   \"output_filename\" " if (@ARGV < 3);
my $dir=$ARGV[0];       # /data/deep_seq9/pacbio/20120720_4xSunny_2xAna_59/A01_1/Analysis_Results/
my $lane=$ARGV[1];      # m120720_155039_42155_c100362222550000001523027210041202
my $fileout=$ARGV[2];
open OUT1,">".$fileout.".subread";
open OUT2,">".$fileout.".ccs";
my $csv1=$dir."/".$lane."_s1_p0.sts.csv";
my $mycsv1=$fileout."_s1_p0.sts.csv";
system("cp $csv1 $mycsv1");

my %Serial;
my %CCS;
open IN, $csv1;
while(<IN>) {
    chomp;
    if (m/Zmw/) {next;}
    my @a=split(",",$_);
    $a[9]=~s/\s//g;
    # marking the ID of successful wells         ReadLength,MedianInsertLength
    if ($a[10] eq 1)  { $Serial{$a[0]}=join("\t",$a[9],$a[11]); }
}
close IN;

# output movie-1
open IN1,$dir."/".$lane."_s1_p0.1.subreads.fasta";
open IN2,$dir."/".$lane."_s1_p0.1.ccs.fasta";
my $id="";
my $seq="";
while(<IN2>) {
    chomp;
    if (m/^>/) {
        if ($seq ne "") {
            my @a=split(" ",$id);
            my @b=split(/\//,$a[0]);
            if (exists $Serial{$b[1]}) {
                $CCS{$b[1]}=1;
                my @c=split("\t",$Serial{$b[1]});
                print OUT2 ">ccs_".$b[1]."_".$c[0]."_0_".length($seq)."_".length($seq),"\n",$seq,"\n";
            }
        }
        #print OUT $id,"\n",$seq,"\n";
        $id=$_;
        $seq="";
    }
    else { $seq=$seq.$_; }
}
if ($seq ne "") {
    my @a=split(" ",$id);
    my @b=split(/\//,$a[0]);
    if (exists $Serial{$b[1]}) {
        $CCS{$b[1]}=1;
        my @c=split("\t",$Serial{$b[1]});
        print OUT2 ">ccs_".$b[1]."_".$c[0]."_0_".length($seq)."_".length($seq),"\n",$seq,"\n";
    }
}
close IN2;
$id="";
$seq="";
while(<IN1>) {
    chomp;
    if (m/^>/) {
        if ($seq ne "") {
            my @a=split(" ",$id);
            my @b=split(/\//,$a[0]);
            if ((exists $Serial{$b[1]}) and (!exists $CCS{$b[1]})) {
                my @c=split("\t",$Serial{$b[1]});
                print OUT1 ">sub_".$b[1]."_".$c[0]."_".$b[2]."_".length($seq),"\n",$seq,"\n";
            }
        }
        #print OUT $id,"\n",$seq,"\n";
        $id=$_;
        $seq="";
    }
    else { $seq=$seq.$_; }
}
if ($seq ne "") {
    my @a=split(" ",$id);
    my @b=split(/\//,$a[0]);
    if ((exists $Serial{$b[1]}) and (!exists $CCS{$b[1]})) {
        my @c=split("\t",$Serial{$b[1]});
        print OUT1 ">sub_".$b[1]."_".$c[0]."_".$b[2]."_".length($seq),"\n",$seq,"\n";
    }
}
close IN1;

# output movie-2
open IN1,$dir."/".$lane."_s1_p0.2.subreads.fasta";
open IN2,$dir."/".$lane."_s1_p0.2.ccs.fasta";
$id="";
$seq="";
while(<IN2>) {
    chomp;
    if (m/^>/) {
        if ($seq ne "") {
            my @a=split(" ",$id);
            my @b=split(/\//,$a[0]);
            if (exists $Serial{$b[1]}) {
                $CCS{$b[1]}=1;
                my @c=split("\t",$Serial{$b[1]});
                print OUT2 ">ccs_".$b[1]."_".$c[0]."_0_".length($seq)."_".length($seq),"\n",$seq,"\n";
            }
        }
        #print OUT $id,"\n",$seq,"\n";
        $id=$_;
        $seq="";
    }
    else { $seq=$seq.$_; }
}
if ($seq ne "") {
    my @a=split(" ",$id);
    my @b=split(/\//,$a[0]);
    if (exists $Serial{$b[1]}) {
        $CCS{$b[1]}=1;
        my @c=split("\t",$Serial{$b[1]});
        print OUT2 ">ccs_".$b[1]."_".$c[0]."_0_".length($seq)."_".length($seq),"\n",$seq,"\n";
    }
}
close IN2;
$id="";
$seq="";
while(<IN1>) {
    chomp;
    if (m/^>/) {
        if ($seq ne "") {
            my @a=split(" ",$id);
            my @b=split(/\//,$a[0]);
            if ((exists $Serial{$b[1]}) and (!exists $CCS{$b[1]})) {
                my @c=split("\t",$Serial{$b[1]});
                print OUT1 ">sub_".$b[1]."_".$c[0]."_".$b[2]."_".length($seq),"\n",$seq,"\n";
            }
        }
        #print OUT $id,"\n",$seq,"\n";
        $id=$_;
        $seq="";
    }
    else { $seq=$seq.$_; }
}
if ($seq ne "") {
    my @a=split(" ",$id);
    my @b=split(/\//,$a[0]);
    if ((exists $Serial{$b[1]}) and (!exists $CCS{$b[1]})) {
        my @c=split("\t",$Serial{$b[1]});
        print OUT1 ">sub_".$b[1]."_".$c[0]."_".$b[2]."_".length($seq),"\n",$seq,"\n";
    }
}
close IN1;

# output movie-3
open IN1,$dir."/".$lane."_s1_p0.3.subreads.fasta";
open IN2,$dir."/".$lane."_s1_p0.3.ccs.fasta";
$id="";
$seq="";
while(<IN2>) {
    chomp;
    if (m/^>/) {
        if ($seq ne "") {
            my @a=split(" ",$id);
            my @b=split(/\//,$a[0]);
            if (exists $Serial{$b[1]}) {
                $CCS{$b[1]}=1;
                my @c=split("\t",$Serial{$b[1]});
                print OUT2 ">ccs_".$b[1]."_".$c[0]."_0_".length($seq)."_".length($seq),"\n",$seq,"\n";
            }
        }
        #print OUT $id,"\n",$seq,"\n";
        $id=$_;
        $seq="";
    }
    else { $seq=$seq.$_; }
}
if ($seq ne "") {
    my @a=split(" ",$id);
    my @b=split(/\//,$a[0]);
    if (exists $Serial{$b[1]}) {
        $CCS{$b[1]}=1;
        my @c=split("\t",$Serial{$b[1]});
        print OUT2 ">ccs_".$b[1]."_".$c[0]."_0_".length($seq)."_".length($seq),"\n",$seq,"\n";
    }
}
close IN2;
$id="";
$seq="";
while(<IN1>) {
    chomp;
    if (m/^>/) {
        if ($seq ne "") {
            my @a=split(" ",$id);
            my @b=split(/\//,$a[0]);
            if ((exists $Serial{$b[1]}) and (!exists $CCS{$b[1]})) {
                my @c=split("\t",$Serial{$b[1]});
                print OUT1 ">sub_".$b[1]."_".$c[0]."_".$b[2]."_".length($seq),"\n",$seq,"\n";
            }
        }
        #print OUT $id,"\n",$seq,"\n";
        $id=$_;
        $seq="";
    }
    else { $seq=$seq.$_; }
}
if ($seq ne "") {
    my @a=split(" ",$id);
    my @b=split(/\//,$a[0]);
    if ((exists $Serial{$b[1]}) and (!exists $CCS{$b[1]})) {
        my @c=split("\t",$Serial{$b[1]});
        print OUT1 ">sub_".$b[1]."_".$c[0]."_".$b[2]."_".length($seq),"\n",$seq,"\n";
    }
}
close IN1;