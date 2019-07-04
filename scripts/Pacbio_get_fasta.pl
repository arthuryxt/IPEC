#!/usr/bin/perl -w
use strict;
die "Usage: $0  \" dir\"  \"serial\/lane number\"   \"output_filename\" " if (@ARGV < 3);
my $dir=$ARGV[0];       # /data/deep_seq9/pacbio/20120720_4xSunny_2xAna_59/A01_1/Analysis_Results/
my $lane=$ARGV[1];      # m120720_155039_42155_c100362222550000001523027210041202
my $fileout=$ARGV[2];
my %Serial;
open OUT,">".$fileout;
open OUT2,">".$fileout.".ccs";
my $csv1=$dir."/".$lane."_s1_p0.sts.csv";
my $mycsv1=$fileout."_s1_p0.sts.csv";
my $csv2=$dir."/".$lane."_s2_p0.sts.csv";
my $mycsv2=$fileout."_s2_p0.sts.csv";
system("cp $csv1 $mycsv1");
system("cp $csv2 $mycsv2");
# output movie-1
open IN, $dir."/".$lane."_s1_p0.sts.csv";
while(<IN>) {
    chomp;
    if (m/Zmw/) {next;}
    my @a=split(",",$_);
    $a[9]=~s/\s//g;
    # marking the ID of successful wells 
    if ($a[9] eq 1)  { $Serial{$a[0]}=1; }
}
close IN;
open IN1,$dir."/".$lane."_s1_p0.fasta";
open IN2,$dir."/".$lane."_s1_p0.ccs.fasta";
my $id="";
my $seq="";
while(<IN1>) {
    chomp;
    if (m/^>/) {
        if ($seq ne "") {
            my @a=split(/\_/,$id);
            my @b=split(/\//,$a[-1]);
            if (exists $Serial{$b[1]}) {
                print OUT ">".$a[-2]."_".$a[-1],"\n",$seq,"\n";
            }
        }
        #print OUT $id,"\n",$seq,"\n";
        $id=$_;
        $seq="";
    }
    else { $seq=$seq.$_; }
}
close IN1;
$id="";
$seq="";
while(<IN2>) {
    chomp;
    if (m/^>/) {
        if ($seq ne "") {
            my @a=split(/\_/,$id);
            my @b=split(/\//,$a[-1]);
            if (exists $Serial{$b[1]}) {
                print OUT2 ">".$a[-2]."_".$b[0]."/".$b[1],"\n",$seq,"\n";
            }
        }
        #print OUT $id,"\n",$seq,"\n";
        $id=$_;
        $seq="";
    }
    else { $seq=$seq.$_; }
}
close IN2;

# output movie-2
my %Serial2;
open IN, $dir."/".$lane."_s2_p0.sts.csv";
while(<IN>) {
    chomp;
    if (m/Zmw/) {next;}
    my @a=split(",",$_);
    $a[9]=~s/\s//g;
    # marking the ID of successful wells 
    if ($a[9] eq 1)  { $Serial2{$a[0]}=1; }
}
close IN;
open IN1,$dir."/".$lane."_s2_p0.fasta";
open IN2,$dir."/".$lane."_s2_p0.ccs.fasta";
$id="";
$seq="";
while(<IN1>) {
    chomp;
    if (m/^>/) {
        if ($id ne "") {
            my @a=split(/\_/,$id);
            my @b=split(/\//,$a[-1]);
            if (exists $Serial2{$b[1]}) {
                print OUT ">".$a[-2]."_".$a[-1],"\n",$seq,"\n";
            }
        }
        $id=$_;
        $seq="";
    }
    else { $seq=$seq.$_; }
}
close IN1;
$id="";
$seq="";
while(<IN2>) {
    chomp;
    if (m/^>/) {
        if ($seq ne "") {
            my @a=split(/\_/,$id);
            my @b=split(/\//,$a[-1]);
            if (exists $Serial{$b[1]}) {
                print OUT2 ">".$a[-2]."_".$b[0]."/".$b[1],"\n",$seq,"\n";
            }
        }
        #print OUT $id,"\n",$seq,"\n";
        $id=$_;
        $seq="";
    }
    else { $seq=$seq.$_; }
}
close IN2