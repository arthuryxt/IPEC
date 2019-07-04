#!/usr/bin/perl -w
use strict;
my $file=$ARGV[0];
open IN,$file;
open OUT,">".$file.".fa";
my $id="";
my $seq="";
while(<IN>) {
    chomp $_;
    if (/^>/) {
        #my @a=split(" ",$_);
        if ($id eq "") {$id=$_;}
        else {
            print OUT $id,"\n";
            print OUT $seq,"\n";
            $id=$_;
            $seq="";
        }
    }
    else {
        if ($seq eq "") {$seq=$_;}
        else{$seq=$seq.$_;}
    }
}
print OUT $id,"\n";
print OUT $seq,"\n";
