#!/usr/bin/perl -w
use strict;
die "Usage: $0   \"input_cluster\"    \"output_basename\"  \"debug defaule=0\" " if (@ARGV < 2);
my $filein1=$ARGV[0];	# A1_rat_20150902_Mn8.polyA.good
my $fileout=$ARGV[1];	# A1_rat_20150902_Mn8.polyA.new
my $debug=0;
if (scalar(@ARGV) > 2) {$debug=$ARGV[2];}
my %uniq;
open IN,$filein1;
open OUT,">".$fileout;
open OUT2,">".$fileout.".loci";
my $GL="";
my $strand="";
my $chr="";
my $f=0;
my %range;

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

while(<IN>) {
	chomp;
	my @a=split("\t",$_);
	if($a[0] eq $GL){
		# check if this transcript overlap any existing range
		$uniq{$GL}{$a[4]}{$a[5]}=join("\t",@a);
		my $fallinside="";
		foreach my $L (sort{$a <=> $b} keys %{$range{$GL}}) {
			my $overlap=check_overlap($a[4],$a[5],$L,$range{$GL}{$L});
			if ($overlap >= 0) {
				if ($fallinside eq "") {$fallinside=$L;}
				else {$fallinside=$fallinside."\t".$L;}
			}
		}
		if ($fallinside ne "") {
			my @b=split("\t",$fallinside);
			my $min=$a[4];
			my $max=$a[5];
			for(my $i=0; $i<scalar(@b); $i++) {
				if ($debug > 0) {
					print $a[4],"\t",$a[5],"\t",$b[$i],"\t",$range{$GL}{$b[$i]},"\n";
				}
				$min=$min < $b[$i] ? $min : $b[$i];
				$max=$max < $range{$GL}{$b[$i]} ? $range{$GL}{$b[$i]} : $max;
				delete $range{$GL}{$b[$i]};
			}
			$range{$GL}{$min}=$max;
		}
		else {
			$range{$GL}{$a[4]}=$a[5];
		}
		if ($debug > 0) { foreach my $L (sort{$a <=> $b} keys %{$range{$GL}}) { print join("\t","save",$L,$range{$GL}{$L}),"\n";} }
	}
	else{
		if($GL ne ""){
			my $Nr=0;
			foreach my $L (keys %{$range{$GL}}) { $Nr++; }
			if ($debug > 0) {print $Nr,"\n";}
			if ($Nr eq 1) {
				my $tcnt=0;
				foreach my $L(sort{$a <=> $b} keys %{$uniq{$GL}}){
					foreach my $R(sort{$a <=> $b}keys %{$uniq{$GL}{$L}}){
						print OUT $uniq{$GL}{$L}{$R},"\n";
						$tcnt++;
					}
				}
				foreach my $L (keys %{$range{$GL}}) {
					print OUT2 join("\t",$chr,$L,$range{$GL}{$L},$GL,$tcnt,$strand),"\n";
				}
			}
			elsif ($Nr > 1) {
				my $tmpNr=0;
				foreach my $L (sort{$a <=> $b} keys %{$range{$GL}}) {
					$tmpNr++;
					my $tcnt=0;
					foreach my $tL(sort{$a <=> $b} keys %{$uniq{$GL}}){
						foreach my $tR(sort{$a <=> $b}keys %{$uniq{$GL}{$tL}}){
							if (($L <= $tL) and ($tR <= $range{$GL}{$L})) {
								my @tmpt=split("\t",$uniq{$GL}{$tL}{$tR});
								$tmpt[0]=$tmpt[0].".".$tmpNr;
								print OUT join("\t",@tmpt),"\n";
								$tcnt++;
							}
						}
					}
					print OUT2 join("\t",$chr,$L,$range{$GL}{$L},$GL.".".$tmpNr,$tcnt,$strand),"\n";
				}
			}
		}
		$f=0;
		$GL=$a[0];
		$chr=$a[2];
		$strand=$a[3];
		$uniq{$GL}{$a[4]}{$a[5]}=join("\t",@a);
		$range{$GL}{$a[4]}=$a[5];
	}
}
{
			my $Nr=0;
			foreach my $L (keys %{$range{$GL}}) { $Nr++; }
			if ($debug > 0) {print $Nr,"\n";}
			if ($Nr eq 1) {
				my $tcnt=0;
				foreach my $L(sort{$a <=> $b} keys %{$uniq{$GL}}){
					foreach my $R(sort{$a <=> $b}keys %{$uniq{$GL}{$L}}){
						print OUT $uniq{$GL}{$L}{$R},"\n";
						$tcnt++;
					}
				}
				foreach my $L (keys %{$range{$GL}}) {
					print OUT2 join("\t",$chr,$L,$range{$GL}{$L},$GL,$tcnt,$strand),"\n";
				}
			}
			elsif ($Nr > 1) {
				my $tmpNr=0;
				foreach my $L (sort{$a <=> $b} keys %{$range{$GL}}) {
					$tmpNr++;
					my $tcnt=0;
					foreach my $tL(sort{$a <=> $b} keys %{$uniq{$GL}}){
						foreach my $tR(sort{$a <=> $b}keys %{$uniq{$GL}{$tL}}){
							if (($L <= $tL) and ($tR <= $range{$GL}{$L})) {
								my @tmpt=split("\t",$uniq{$GL}{$tL}{$tR});
								$tmpt[0]=$tmpt[0].".".$tmpNr;
								print OUT join("\t",@tmpt),"\n";
								$tcnt++;
							}
						}
					}
					print OUT2 join("\t",$chr,$L,$range{$GL}{$L},$GL.".".$tmpNr,$tcnt,$strand),"\n";
				}
			}
}


