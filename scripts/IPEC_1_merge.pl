#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"out_file_basename\"    \"mod_file_list\"    " if (@ARGV ne 2);
# merge mod file and report final corrected fa
# run after : perl ~/commoncode/IPEC_1.pl ./1_index/Spret_LD4.101.ccs ./2_mapping/LD4.101.ccs.SFL_000.bt2aln ./3_LD4/art_LD4.101.ccs.SFL_000.fa 2

my @time1 = localtime(time);
my @time2;

my $fileout=$ARGV[0];
my $filelist=$ARGV[1];

my %Uniq;
open INL, $filelist;
my $filecnt=0;
while(<INL>) {
	chomp;
	my $filein=$_;
	print "reading mod file  : $filein\n";
	@time2 = localtime(time);
	print join("\t"," running time: ",$time2[2]-$time1[2],"h",$time2[1]-$time1[1],"min",$time2[0]-$time1[0],"sec"),"\n";
    $filecnt++;
	open IN, $filein;
	while(<IN>) {
		chomp;
		my @a=split("\t",$_);
		my @b=split("\;",$a[5]);
		for(@b) {
			my @c=split("\=",$_);
			if (($filecnt > 1) and ($c[1] > int($c[1]))) {
				$Uniq{$a[0]}{$a[1]}{$c[0]}+=($c[1]-1.5);
		    }
			else { $Uniq{$a[0]}{$a[1]}{$c[0]}+=$c[1]; }
		}
	}
	close IN;
}
print "reporting...\n";

open OUT,">".$fileout.".fa";
open OUT1,">".$fileout.".cov";
open OUT11,">".$fileout.".cov1";
open OUT2,">".$fileout.".mod"; 
foreach my $id (keys %Uniq) {
    my $cov="";
    my $seq="";
    my $ccov=0;
    foreach my $pos (sort{$a <=> $b} keys %{$Uniq{$id}}) {
	my $max=0;
	my $maxNT="";
	my $info="";
	my $sum=0;
	foreach my $NT (sort{$Uniq{$id}{$pos}{$b} <=> $Uniq{$id}{$pos}{$a}} keys %{$Uniq{$id}{$pos}}) {
	    if ($maxNT eq "") {
		$max=$Uniq{$id}{$pos}{$NT};
		$maxNT=$NT;
		$sum=$max;
		$info=$NT."=".$max;
	    }
	    else {
		$info=$info.";".$NT."=".$Uniq{$id}{$pos}{$NT};
		$sum+=$Uniq{$id}{$pos}{$NT};
	    }
	}
	print OUT2 join("\t",$id,$pos,$max,$sum,sprintf("%.3f",$max/$sum),$info),"\n"; 
	if (($maxNT ne "Del") and ($maxNT ne "_") and ($maxNT ne "")) {
	    $seq=$seq.$maxNT;
	    my $tmp_Nr=length($maxNT);
    	if ($sum eq 1.5) {
		    for(my $tmp=0; $tmp<$tmp_Nr; $tmp++) { $cov=$cov."0"; }
	    }
	    else {
		    for(my $tmp=0; $tmp<$tmp_Nr; $tmp++) { $cov=$cov."1"; $ccov++; }
	    }
    }
	#if ($maxNT ne "Del")  { $seq=$seq.$maxNT;  }
    }
    print OUT ">".$id,"\n",$seq,"\n";
    print OUT1 ">".$id,"\n",$cov,"\n";
    print OUT11 $id,"\t",$ccov,"\t",length($cov),"\t",sprintf("%.2f",$ccov/length($cov)),"\n";
}
close OUT;
close OUT1;

@time2 = localtime(time);
print join("\t"," running time: ",$time2[2]-$time1[2],"h",$time2[1]-$time1[1],"min",$time2[0]-$time1[0],"sec"),"\n";
