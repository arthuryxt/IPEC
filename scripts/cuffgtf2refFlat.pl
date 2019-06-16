#!/usr/bin/perl -w
use strict;
# 2015-01-19 by Arthur
die "Usage: $0   \"cuff_gtf\"   \"output_refFlat\"   \(DEBUG==1\)" if (@ARGV < 2);
my $filein=$ARGV[0];
my $fileout=$ARGV[1];
my $debug=0;
if (scalar(@ARGV) > 2) { $debug=$ARGV[2]; }

open IN, $filein;
open OUT,">".$fileout;

my $locus="";
my $transcrpt="";
my $gene_name="";
my $chr="";
my $strand="";
my $exonNr=0;
my @Start="";
my @End="";
while (<IN>) {
    chomp;
    if (m/^#/) { next; }
    my @a=split("\t",$_);
    $a[8]=~s/\s//g;
    $a[8]=~s/\;//g;
    my @b=split(/\"/,$a[8]);
    my $tlocus="";
    my $ttranscrpt="";
    my $tgene_name="";
    my $tchr=$a[0];
    if ($tchr!~m/chr/) { $tchr="chr".$tchr; }
    my $tstrand=$a[6];
    my $tstart=$a[3]-1;
    my $tend=$a[4];
    my $nr=scalar(@b);
    for(my $i=0; $i<$nr; $i=$i+2) {
        if ($b[$i] eq "gene_id") { $tlocus=$b[$i+1]; }
        if ($b[$i] eq "transcript_id") { $ttranscrpt=$b[$i+1]; }
        if ($b[$i] eq "gene_name") { $tgene_name=$b[$i+1]; }
    }
    if ($transcrpt eq "") {
        $locus=$tlocus;
        $transcrpt=$ttranscrpt;
        $gene_name=$tgene_name;
        $chr=$tchr;
        $strand=$tstrand;
        $Start[$exonNr]=$tstart;
        $End[$exonNr]=$tend;
        $exonNr++;
    }
    elsif ($transcrpt eq $ttranscrpt) {
        $Start[$exonNr]=$tstart;
        $End[$exonNr]=$tend;
        $exonNr++;
    }
    else {
        my $starts=$Start[0];
        for(my $i=1; $i<$exonNr; $i++) { $starts=$starts.",".$Start[$i]; }
        my $ends=$End[0];
        for(my $i=1; $i<$exonNr; $i++) { $ends=$ends.",".$End[$i]; }
        print OUT join("\t",$locus,$transcrpt,$chr,$strand,$Start[0],$End[$exonNr-1],$Start[0],$Start[0],$exonNr,$starts,$ends,$gene_name),"\n";
        
        $locus=$tlocus;
        $transcrpt=$ttranscrpt;
        $gene_name=$tgene_name;
        $chr=$tchr;
        $strand=$tstrand;
        $exonNr=0;
        $Start[$exonNr]=$tstart;
        $End[$exonNr]=$tend;
        $exonNr++;
    }
}
