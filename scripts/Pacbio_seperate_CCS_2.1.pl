#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"subreads.t15\"  \"ReadOfInsert\"  \"output_subreads\"  \"output_ccs\"  \"\(optional\)label_prefix_header\"" if (@ARGV < 4);
my $filein1=$ARGV[0];
my $filein2=$ARGV[1];
my $fileout1=$ARGV[2];
my $fileout2=$ARGV[3];
my $label=0;
if (scalar(@ARGV) > 4){ $label=$ARGV[4];}

my %SUB;
open(IN, $filein1) or die "cannot locate input_subread.fa";
while (<IN>) {
	chomp;
	s/^>//;
	my @a=split("_",$_);
	my $seq=<IN>;
	chomp $seq;
	#				   0          1          2          3          4          5          6
	$SUB{$a[1]}{$a[3]}{$a[0]."\t".$a[2]."\t".$a[4]."\t".$a[5]."\t".$a[6]."\t".$a[7]."\t".$a[8]}=$seq;
}
close IN;

my %CCS;
open(IN2, $filein2) or die "cannot locate ReadOfInsert";
close IN2;
my $command="perl /home/arthur/commoncode/flat454_fullheader_no_length.pl $filein2";
system($command);
my $filein2fa=$filein2.".fa";
open IN2,$filein2fa;

open(OUT2, ">".$fileout2) or die "cannot open output_ccs";
while (<IN2>) {
	chomp;
	s/^>//;
	my @a=split("\/",$_);
	my $seq=<IN2>;
	chomp $seq;
	my $len=length($seq);
	$CCS{$a[1]}=1;
	if ($label eq 0) {
		print OUT2 ">ccs_".$a[1]."_".$len."_0_".$len."_".$len,"\n",$seq,"\n";
	}
	else {
		print OUT2 ">".$label."_ccs_".$a[1]."_".$len."_0_".$len."_".$len,"\n",$seq,"\n";
	}
}
close IN2;
close OUT2;
$command="rm -f $filein2fa";
system($command);

open(OUT1, ">".$fileout1) or die "cannot open output_subreads";
foreach my $id (sort{$a <=> $b} keys %SUB) {
	if (exists $CCS{$id}) {
		next;
	}
	foreach my $start (sort{$a <=> $b} keys %{$SUB{$id}}) {
		foreach my $string (sort keys %{$SUB{$id}{$start}}) {
			my @b=split("\t",$string);
			if ($label eq 0) {
				print OUT1 ">".$b[0]."_".$id."_".$b[1]."_".$start."_".$b[2]."_".$b[3]."_".$b[4]."_".$b[5]."_".$b[6],"\n",$SUB{$id}{$start}{$string},"\n";
			}
			else {
				print OUT1 ">".$label."_".$b[0]."_".$id."_".$b[1]."_".$start."_".$b[2]."_".$b[3]."_".$b[4]."_".$b[5]."_".$b[6],"\n",$SUB{$id}{$start}{$string},"\n";
			}
		}
	}
}
close OUT1;