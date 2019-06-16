#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"reference.fa\"    \"HQ.alignment\"    \"output\"   \"\(optional\)debug\" " if (@ARGV < 3);
# run after : qsub -N PBM -l h_vmem=4G -pe smp 6 -cwd -V -b y '~/bin/bowtie2-2.0.2/bowtie2 -p 20 -q --phred64 --local -D 40, -R 3, -N 1 -L 20 -i S,1,0.50 -k 10 --no-head -x EC_step1.subreads EC_step1.illumina.fq.uniq > EC_step1.illumina.bt2aln'

my @time1 = localtime(time);
my @time2;

my $filein1=$ARGV[0];	# EC_step1.subreads
my $filein2=$ARGV[1];	# EC_step1.illumina.bt2aln
my $fileout=$ARGV[2];	# EC_step1.IPEC_1
my $debug=0;
if (scalar(@ARGV) > 3) {$debug=$ARGV[3];}
my %SEQ;
my %Uniq;
open IN,$filein1;
while(<IN>) {
    chomp;
    s/^>//;
    my $id=$_;
    my $seq=<IN>;
    chomp $seq;
    $SEQ{$id}=$seq;
    my @a=split(//,$seq);
    my $Nr=length($seq);
    for(my $i=0; $i<$Nr; $i++) { $Uniq{$id}{$i}{$a[$i]}=1.5; }
}
close IN;

sub magic($$$) {
    my $cigar=shift;
    my $md=shift;
    my $seq=shift;
    # parse $cigar
    my @cigar_move=($cigar=~m/\d+(\w)/g);
    my @cigar_value=($cigar=~m/(\d+)\w/g);
    my $cigar_Nr=scalar(@cigar_move);
    my @cigar_master;
    my @cigar_seq;
    my $cnt=0;
    my $pos_on_seq=0;
    for(my $i=0; $i<$cigar_Nr; $i++) {
	for(my $j=0; $j<$cigar_value[$i]; $j++) {
	    $cigar_master[$cnt]=$cigar_move[$i];
            if ($cigar_move[$i] eq "D") {$cigar_seq[$cnt]="_";}
    	    else {$cigar_seq[$cnt]=substr($seq,$pos_on_seq,1); $pos_on_seq++;}
	    $cnt++;
	}
    }
    # parse md
    $md=~m/(\d*)/;
    my $md_0=$1;
    my @md_array=($md=~m/([A-Z]|\^[A-Z]+)(\d+)/g);
    my @md_master=@cigar_master;
    my $count=0;
    my $tmp=$md_0;
    while($tmp > 0) {
	if ($md_master[$count] eq "M") {$tmp--;}
	$count++; 
    }
    for(@md_array) {
	my $move=$_;
	if($move=~m/^\^/) {
	    # insertion on the reference, already marked as "D" in the master_array
	    $move=~s/\^//;
    	    $count+=length($move);
	}
	elsif ($move=~m/^\d+$/) {
	    # match and insertion
	    my $tmp=$move;
	    while(($tmp > 0) and ($count < $cnt)) {
	        if (($md_master[$count] eq "M") or ($md_master[$count] eq "D")) { $tmp--; }
    	        $count++; 
	    }
	}
	else {
	    # mismatch
	    if ($cigar_seq[$count] ne "_") { $md_master[$count]=$cigar_seq[$count]; }
	    $count++;
	}
    }
    return join("",@cigar_seq)."\t".join("",@md_master);
}

sub update($$$$$) {
    my $refid=shift;
    my $weight=shift;
    my $start=shift;
    my @cigar_seq=split//, shift;
    my @master=split//, shift;
    my $cnt=scalar(@master);
    my $pos_on_ref=$start;
    if ($debug eq 998) { print $refid,"\t",$weight,"\t",$start,"\n","CIGAR  ",join("",@cigar_seq),"\n","master ",join("",@master),"\n";} 
    for(my $i=0; $i<$cnt; $i++) {
	# clean master_array, will only start with "M" out of [MIDATCG]
	if ((($i+1)<$cnt) and ($master[$i+1] ne "I")) {
	    if ($master[$i] eq "D") { $Uniq{$refid}{$pos_on_ref}{"Del"}+=$weight; $pos_on_ref++;}
	    elsif ($master[$i] eq "M") {
		$Uniq{$refid}{$pos_on_ref}{$cigar_seq[$i]}+=$weight; $pos_on_ref++;
		if ($debug eq 1) { if ($cigar_seq[$i]=~m/\_/) {print $refid,"\t",$weight,"\t",$start,"\n","CIGAR  ",join("",@cigar_seq),"\n","master ",join("",@master),"\n";} }
	    }
	    else {$Uniq{$refid}{$pos_on_ref}{$master[$i]}+=$weight; $pos_on_ref++;}
	}
	else {
	    my $NT=$cigar_seq[$i];
	    if ($cigar_seq[$i] eq "_") {$NT="";}
	    if ($master[$i]=~m/[ATCG]/) {$NT=$master[$i];}
	    while((($i+1)<$cnt) and ($master[$i+1] eq "I")){
		if ($cigar_seq[$i+1] ne "_") { $NT=$NT.$cigar_seq[$i+1]; }
		#$NT=$NT.$cigar_seq[$i+1];
		$i++;
	    }
	    if ($debug eq 1) { if ($NT=~m/\_/) {print $refid,"\t",$weight,"\t",$start,"\n","CIGAR  ",join("",@cigar_seq),"\n","master ",join("",@master),"\n";} }
	    $Uniq{$refid}{$pos_on_ref}{$NT}+=$weight;
	    $pos_on_ref++;
	}
    }
}

open IN,$filein2;
#open OUTp,">".$fileout.".tmp";
while(<IN>) {
    chomp;
    if ((m/^@/) or (m/^#/)) {next;}
    my @a=split("\t",$_);
    #	0					1	2				3	4	5					6	7	8	9	10	11...
    #	HWI-ST540.153.2.1211.10066.96626/1__1	256	s1_p0/100_7028_0_3041_3041	43	11	10M1I9M1D13M1D6M1D16M2D35M1I7M1I22M	*	0	0	seq	qua	AS:i:96	XS:i:138	MD:Z:3A15^C10A2^G6^C7G0T0A6^AA1A0C4C11C11A0A0C2A1C0G2A1T1A12G2A1
    if ($a[1] eq 4) { next; }
    if (!exists $SEQ{$a[2]}) {next;}
    my $refid=$a[2];
    my $FLAG=$a[1];
    my $start=$a[3]-1;
    my $CIGAR=$a[5];
    my $sequence=$a[9];
    my $MD="";
    my $Nr=scalar(@a);
    for(my $i=11; $i<$Nr; $i++) {
	if($a[$i]=~m/MD/) { my @tmp=split(/\:/,$a[$i]); $MD=$tmp[2]; last; }
    }
    my $count=1;
    if ($a[0]=~m/\_\_/){ my @tmp=split(/\_\_/,$a[0]); $count=$tmp[-1]; }
    my @info=split("\t",magic($CIGAR,$MD,$sequence));
    # quality control for soft-clip
    $info[1]=~m/(S*)([^S]*)(S*)/;
    my $left=length($1);
    my $mid=length($2);
    my $right=length($3);
    #print OUTp join("\t",$refid,$left,$mid,$right,$start),"\n",$info[0],"\n",$info[1],"\n";
    if (($left+$right) > $mid) { next; }
    my $str1=substr($info[0],$left,$mid);
    my $str2=substr($info[1],$left,$mid);
    #print OUTp join("\n",$refid,$count."\t".$a[0],$start,$str1,$str2),"\n";
    update($refid,$count,$start,$str1,$str2);
}
close IN;

open OUT,">".$fileout;
open OUT1,">".$fileout.".cov";
open OUT11,">".$fileout.".cov1";
if ($debug eq 2) { open OUT2,">".$fileout.".mod"; }
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
	if ($debug eq 2) { print OUT2 join("\t",$id,$pos,$max,$sum,sprintf("%.3f",$max/$sum),$info),"\n"; }
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
if ($debug eq 2) { close OUT2; }

@time2 = localtime(time);
print join("\t"," running time: ",$time2[2]-$time1[2],"h",$time2[1]-$time1[1],"min",$time2[0]-$time1[0],"sec"),"\n";
