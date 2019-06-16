#!/usr/bin/perl -w
use strict;
die "Usage: $0   \"input_border\"    \"input_wiggle\"    \"output_basename\"   \"window_size==1000\"    \"cutoff==0.05\"   \"mincov==2\"   \(strandness==+\)  \(debug==0\)" if (@ARGV < 3);
my $filein1=$ARGV[0];		# A1_rat_20150902_Mn8_fix.overlap.RefSeq.3.border
my $filein2=$ARGV[1];		# r04.rmd.neur1_minus.wig
my $fileout=$ARGV[2];		# r04.rmd.neur1_minus.checkFLT
my $window=1000;
if (scalar(@ARGV) > 3) { $window=$ARGV[3]; }
my $minpctg=0.05;
if (scalar(@ARGV) > 4) { $minpctg=$ARGV[4]; }
my $mincov=2;
if (scalar(@ARGV) > 5) { $mincov=$ARGV[5]; }
my $strand="+";
if (scalar(@ARGV) > 6) { $strand=$ARGV[6]; }
die "strandness should either be + or -\n" if (($strand ne "+") and ($strand ne "-"));
my $debug=0;
if (scalar(@ARGV) > 7) { $debug=$ARGV[7]; }
my $binsize=10000;
my %uniq;
my %GL;
my %GLcov;
open IN1,$filein1;
my %uniq0;
while(<IN1>) {
	chomp;
	my @a=split("\t",$_);
	# chr1    3369291 3428904 GL.8693 +       3369291 3428904
	$a[0]=~s/chr//;
	$uniq0{$a[0]}{$a[4]}{$a[5]}=join("\t",@a);
	#my $bin1=int(($a[5]-$window)/$binsize);
	#my $bin2=int(($a[6]+$window)/$binsize);
	#for(my $i=$bin1; $i<=$bin2; $i++){
	#	$uniq{$a[0]}{$a[4]}{$i}{$a[5]-$window}=join("\t",@a);
	#}
	#$GL{$a[3]}=join("\t",$_);
}
close IN1;
foreach my $chr (sort{$a cmp $b} keys %uniq0) {
	foreach my $strand (sort{$a cmp $b} keys %{$uniq0{$chr}}) {
		my $cnt=0;
		my @pos;
		foreach my $start (sort{$a <=> $b} keys %{$uniq0{$chr}{$strand}}) {
			$pos[$cnt]=$start;
			$cnt++;
		}
		# process the first
		{
			my @tmp0=split("\t",$uniq0{$chr}{$strand}{$pos[0]});
			my @tmp1=split("\t",$uniq0{$chr}{$strand}{$pos[1]});
			my $dist=$tmp1[5] - $tmp0[6];
			if ($dist > 2*$window) {
				my $info=join("\t",$tmp0[0],$tmp0[1],$tmp0[2],$tmp0[3],$tmp0[4],$tmp0[5]-$window,$tmp0[6]+$window,$window);
				$GL{$tmp0[3]}=$info;
				my $bin1=int(($tmp0[5]-$window)/$binsize);
				my $bin2=int(($tmp0[6]+$window)/$binsize);
				for(my $i=$bin1; $i<=$bin2; $i++){ $uniq{$chr}{$strand}{$i}{$tmp0[5]-$window}=$info; }
			}
			else {
				$dist=int($dist/3);
				my $info=join("\t",$tmp0[0],$tmp0[1],$tmp0[2],$tmp0[3],$tmp0[4],$tmp0[5]-$dist,$tmp0[6]+$dist,$dist);
				$GL{$tmp0[3]}=$info;
				my $bin1=int(($tmp0[5]-$dist)/$binsize);
				my $bin2=int(($tmp0[6]+$dist)/$binsize);
				for(my $i=$bin1; $i<=$bin2; $i++){ $uniq{$chr}{$strand}{$i}{$tmp0[5]-$dist}=$info; }
			}
		}
		# process the middle ones
		for(my $i=1; $i<($cnt-1); $i++) {
			my @tmp0=split("\t",$uniq0{$chr}{$strand}{$pos[$i-1]});
			my @tmp1=split("\t",$uniq0{$chr}{$strand}{$pos[$i]});
			my @tmp2=split("\t",$uniq0{$chr}{$strand}{$pos[$i+1]});
			my $dist1=$tmp1[5] - $tmp0[6];
			my $dist2=$tmp2[5] - $tmp1[6];
			my $dist=$dist1 < $dist2 ? $dist1 : $dist2;
			if ($dist > 2*$window) {
				my $info=join("\t",$tmp1[0],$tmp1[1],$tmp1[2],$tmp1[3],$tmp1[4],$tmp1[5]-$window,$tmp1[6]+$window,$window);
				$GL{$tmp1[3]}=$info;
				my $bin1=int(($tmp1[5]-$window)/$binsize);
				my $bin2=int(($tmp1[6]+$window)/$binsize);
				for(my $i=$bin1; $i<=$bin2; $i++){ $uniq{$chr}{$strand}{$i}{$tmp1[5]-$window}=$info; }
			}
			else {
				$dist=int($dist/3);
				my $info=join("\t",$tmp1[0],$tmp1[1],$tmp1[2],$tmp1[3],$tmp1[4],$tmp1[5]-$dist,$tmp1[6]+$dist,$dist);
				$GL{$tmp1[3]}=$info;
				my $bin1=int(($tmp1[5]-$window)/$binsize);
				my $bin2=int(($tmp1[6]+$window)/$binsize);
				for(my $i=$bin1; $i<=$bin2; $i++){ $uniq{$chr}{$strand}{$i}{$tmp1[5]-$dist}=$info; }
			}
		}
		# process the last
		{
			my @tmp0=split("\t",$uniq0{$chr}{$strand}{$pos[$cnt-2]});
			my @tmp1=split("\t",$uniq0{$chr}{$strand}{$pos[$cnt-1]});
			my $dist=$tmp1[5] - $tmp0[6];
			if ($dist > 2*$window) {
				my $info=join("\t",$tmp1[0],$tmp1[1],$tmp1[2],$tmp1[3],$tmp1[4],$tmp1[5]-$window,$tmp1[6]+$window,$window);
				$GL{$tmp1[3]}=$info;
				my $bin1=int(($tmp1[5]-$window)/$binsize);
				my $bin2=int(($tmp1[6]+$window)/$binsize);
				for(my $i=$bin1; $i<=$bin2; $i++){ $uniq{$chr}{$strand}{$i}{$tmp1[5]-$window}=$info; }
			}
			else {
				$dist=int($dist/3);
				my $info=join("\t",$tmp1[0],$tmp1[1],$tmp1[2],$tmp1[3],$tmp1[4],$tmp1[5]-$dist,$tmp1[6]+$dist,$dist);
				$GL{$tmp1[3]}=$info;
				my $bin1=int(($tmp1[5]-$window)/$binsize);
				my $bin2=int(($tmp1[6]+$window)/$binsize);
				for(my $i=$bin1; $i<=$bin2; $i++){ $uniq{$chr}{$strand}{$i}{$tmp1[5]-$dist}=$info; }
			}
		}
	}
}
my %MaxCovInPB;
my %MaxCovOutPB;
my %CntInPB;
my %CntOutPB;
sub check_overlap($$$$){
	if (($_[1] <= $_[2]) or ($_[3] <= $_[0])) { return(-1); }		# 01 23, 23 01
	if ($_[0] < $_[2]) {	# 0 2 1 3
		if ($_[3] < $_[1]) { return(abs($_[3] - $_[2] + 1));}	# 0 2 3 1
		else{ return(abs($_[1] - $_[2]) + 1);}					# 0 2 1 3
	}
	else {
		if ($_[1] < $_[3]) { return(abs($_[1] - $_[0] + 1));}	# 2 0 1 3 
		else{ return(abs($_[3] - $_[0]) + 1);}					# 2 0 3 1
	}
}
open IN2,$filein2;
while(<IN2>) {
	chomp;
	if (m/^track/) {next;}
	if (m/^#/) {next;}
	my @a=split("\t",$_);
	$a[0]=~s/chr//;
	# chr11   185184  185336  1
	my $bin=int($a[1]/$binsize);
	my $f=0;
	if (exists $uniq{$a[0]}{$strand}{$bin}) {
		foreach my $start (keys %{$uniq{$a[0]}{$strand}{$bin}}) {
			my @b=split("\t",$uniq{$a[0]}{$strand}{$bin}{$start});
			my $overlap1=check_overlap($a[1],$a[2],$b[1],$b[2]);
			my $overlap2=check_overlap($a[1],$a[2],$b[5],$b[6]);
			if ($debug > 0) {print join("\t",@b),"\t",join("\t",@a),"\t",$overlap1,"\t",$overlap2,"\n";}
			if (($overlap1 eq -1) and ($overlap2 eq -1)) {}
			elsif(($overlap1 eq -1) and ($overlap2 > 0)) {
				# outside of PB 
				if (exists $MaxCovOutPB{$b[3]}) {
					if ($MaxCovOutPB{$b[3]} < $a[3]) { $MaxCovOutPB{$b[3]}=$a[3]; }
				}
				else {
					$MaxCovOutPB{$b[3]}=$a[3];
				}
				$CntOutPB{$b[3]}+=($a[2] - $a[1]) * $overlap2;
				$f=1;
				for(my $i=$a[1]; $i<$a[2]; $i++){ $GLcov{$b[3]}{$i}=$a[3];}
			}
			elsif($overlap1 > 0) {
				# inside of PB 
				if (exists $MaxCovInPB{$b[3]}) {
					if ($MaxCovInPB{$b[3]} < $a[3]) { $MaxCovInPB{$b[3]}=$a[3]; }
					#if ($MaxCovOutPB{$b[3]} < $a[3]) { $MaxCovOutPB{$b[3]}=$a[3]; }
				}
				else {
					$MaxCovInPB{$b[3]}=$a[3];
					#$MaxCovOutPB{$b[3]}=$a[3];
				}
				$CntInPB{$b[3]}+=($a[2] - $a[1]) * $overlap1;
				#$CntOutPB{$b[3]}+=($a[2] - $a[1]) * $overlap1;
				$f=1;
				for(my $i=$a[1]; $i<$a[2]; $i++){ $GLcov{$b[3]}{$i}=$a[3];}
			}
		}
	}
	if (($f eq 0) and ($uniq{$a[0]}{$strand}{$bin-1})) {
		foreach my $start (keys %{$uniq{$a[0]}{$strand}{$bin-1}}) {
			my @b=split("\t",$uniq{$a[0]}{$strand}{$bin-1}{$start});
			my $overlap1=check_overlap($a[1],$a[2],$b[1],$b[2]);
			my $overlap2=check_overlap($a[1],$a[2],$b[5],$b[6]);
			if ($debug > 0) {print join("\t",@b),"\t",join("\t",@a),"\t",$overlap1,"\t",$overlap2,"\n";}
			if (($overlap1 eq -1) and ($overlap2 eq -1)) {}
			elsif(($overlap1 eq -1) and ($overlap2 > 0)) {
				# outside of PB 
				if (exists $MaxCovOutPB{$b[3]}) {
					if ($MaxCovOutPB{$b[3]} < $a[3]) { $MaxCovOutPB{$b[3]}=$a[3]; }
				}
				else {
					$MaxCovOutPB{$b[3]}=$a[3];
				}
				$CntOutPB{$b[3]}+=($a[2] - $a[1]) * $overlap2;
				$f=1;
				for(my $i=$a[1]; $i<$a[2]; $i++){ $GLcov{$b[3]}{$i}=$a[3];}
			}
			elsif($overlap1 > 0) {
				# inside of PB 
				if (exists $MaxCovInPB{$b[3]}) {
					if ($MaxCovInPB{$b[3]} < $a[3]) { $MaxCovInPB{$b[3]}=$a[3]; }
					#if ($MaxCovOutPB{$b[3]} < $a[3]) { $MaxCovOutPB{$b[3]}=$a[3]; }
				}
				else {
					$MaxCovInPB{$b[3]}=$a[3];
					#$MaxCovOutPB{$b[3]}=$a[3];
				}
				$CntInPB{$b[3]}+=($a[2] - $a[1]) * $overlap1;
				#$CntOutPB{$b[3]}+=($a[2] - $a[1]) * $overlap1;
				$f=1;
				for(my $i=$a[1]; $i<$a[2]; $i++){ $GLcov{$b[3]}{$i}=$a[3];}
			}
		}
	}
	if (($f eq 0) and ($uniq{$a[0]}{$strand}{$bin+1})) {
		foreach my $start (keys %{$uniq{$a[0]}{$strand}{$bin+1}}) {
			my @b=split("\t",$uniq{$a[0]}{$strand}{$bin+1}{$start});
			my $overlap1=check_overlap($a[1],$a[2],$b[1],$b[2]);
			my $overlap2=check_overlap($a[1],$a[2],$b[5],$b[6]);
			if ($debug > 0) {print join("\t",@b),"\t",join("\t",@a),"\t",$overlap1,"\t",$overlap2,"\n";}
			if (($overlap1 eq -1) and ($overlap2 eq -1)) {}
			elsif(($overlap1 eq -1) and ($overlap2 > 0)) {
				# outside of PB 
				if (exists $MaxCovOutPB{$b[3]}) {
					if ($MaxCovOutPB{$b[3]} < $a[3]) { $MaxCovOutPB{$b[3]}=$a[3]; }
				}
				else {
					$MaxCovOutPB{$b[3]}=$a[3];
				}
				$CntOutPB{$b[3]}+=($a[2] - $a[1]) * $overlap2;
				$f=1;
				for(my $i=$a[1]; $i<$a[2]; $i++){ $GLcov{$b[3]}{$i}=$a[3];}
			}
			elsif($overlap1 > 0) {
				# inside of PB 
				if (exists $MaxCovInPB{$b[3]}) {
					if ($MaxCovInPB{$b[3]} < $a[3]) { $MaxCovInPB{$b[3]}=$a[3]; }
					#if ($MaxCovOutPB{$b[3]} < $a[3]) { $MaxCovOutPB{$b[3]}=$a[3]; }
				}
				else {
					$MaxCovInPB{$b[3]}=$a[3];
					#$MaxCovOutPB{$b[3]}=$a[3];
				}
				$CntInPB{$b[3]}+=($a[2] - $a[1]) * $overlap1;
				#$CntOutPB{$b[3]}+=($a[2] - $a[1]) * $overlap1;
				$f=1;
				for(my $i=$a[1]; $i<$a[2]; $i++){ $GLcov{$b[3]}{$i}=$a[3];}
			}
		}
	}
}
close IN2;
open OUT,">".$fileout;
foreach my $id (keys %GL) {
	if ((exists $MaxCovInPB{$id}) or ((exists $MaxCovOutPB{$id}))) {
		my $CovIn=1;
		if (exists $MaxCovInPB{$id}) { $CovIn=$MaxCovInPB{$id}; }
		my $CovOut=0;
		if (exists $MaxCovOutPB{$id}) { $CovOut=$MaxCovOutPB{$id}; }
		my $CntIn=1;
		if (exists $CntInPB{$id}) { $CntIn=$CntInPB{$id}; }
		my $CntOut=0;
		if (exists $CntOutPB{$id}) { $CntOut=$CntOutPB{$id}; }
		my $maxcov=$CovIn > $CovOut ? $CovIn : $CovOut;
		my $newLeft=-1;
		foreach my $pos (sort{$a <=> $b} keys %{$GLcov{$id}}) {
			if ($GLcov{$id}{$pos} > $minpctg * $maxcov) { $newLeft=$pos; last; }
		}
		my $newRight=-1;
		foreach my $pos (sort{$b <=> $a} keys %{$GLcov{$id}}) {
			if ($GLcov{$id}{$pos} > $minpctg * $maxcov) { $newRight=$pos; last; }
		}
		print OUT "chr".$GL{$id},"\t",$CovIn,"\t",$CovOut,"\t",sprintf("%.3f",$CovOut/$CovIn),"\t",$CntIn,"\t",$CntOut,"\t",sprintf("%.3f",$CntOut/$CntIn),"\t",$newLeft,"\t",$newRight,"\n";
	}
}
