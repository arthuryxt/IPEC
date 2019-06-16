#!/usr/bin/perl -w
use strict;
# 2015-01-13 by Arthur
die "Usage: $0   \"input_CAGE_bed6\"    \"input_refFlat\"    \"output\"    \"window\(default==50nt\)\"  \"debug > 0\"  " if (@ARGV < 3);
my $filein1=$ARGV[0];
my $filein2=$ARGV[1];
my $fileout=$ARGV[2];
my $window=50;
if (scalar(@ARGV) > 3) {$window=$ARGV[3];}
my $debug=0;
if (scalar(@ARGV) > 4) {$debug=$ARGV[4];}


my %uniq;
open IN1,$filein1;
while(<IN1>) {
	chomp;
	if (m/^@/) {next;}
	if (m/^#/) {next;}
	if (m/^track/) {next;}
	my @a=split("\t",$_);
	# chr10   28312   28325   TC      44      +
	if ($a[5] eq "+") {
		my $bin1=int($a[1]/$window);
		my $bin2=int($a[2]/$window);
		for(my $i=$bin1; $i<=($bin2+1); $i++) { $uniq{$a[0]}{$a[5]}{$i}{$a[1]}=join("\t",$a[1],$a[2],$a[4]); }
	}
	else {
		my $bin1=int($a[1]/$window);
		my $bin2=int($a[2]/$window);
		for(my $i=$bin1; $i<=($bin2+1); $i++) { $uniq{$a[0]}{$a[5]}{$i}{$a[2]}=join("\t",$a[1],$a[2],$a[4]); }
	}
}
close IN1;

open IN2, $filein2;
open OUT,">".$fileout;
open OUT1,">".$fileout.".dist";	# distance of (TSS - CAGE) on RNA's direction
while(<IN2>) {
	chomp;
	my @a=split("\t",$_);
	#if ($a[2]!~m/^chr/) {$a[2]="chr".$a[2];}
	# GL.1    RSN_LD4.110_ccs_139769_4019_0_4019_4019_SUB0_pNAp_28    chr10   +       10039503        10045902        10039503        10039503  2       10039503,10042783       10039667,10045902       RSN_LD4.110_ccs_139769_4019_0_4019_4019_SUB0_pNAp_28
	if ($a[8] >= 1) {
		if ($a[3] eq "+") {
			my $bin=int($a[4]/$window);
			if ($debug > 0) { print "bin = ",$bin,"\n"; }
			my $flag=0;
			my $info="";
			if (exists $uniq{$a[2]}{$a[3]}{$bin}) {
				foreach my $pos(keys %{$uniq{$a[2]}{$a[3]}{$bin}}) {
					my @b=split("\t",$uniq{$a[2]}{$a[3]}{$bin}{$pos});
					my $distL=$a[4]-$b[0];
					my $distR=$a[4]-$b[1];
					if ($debug > 0) { print join("\t",$a[1],$a[2],$a[3],$a[4],$a[5],"bin",$bin,$distL,$distR),"\n"; }
					if ($distL * $distR <= 0) {
						# 5'end sit inside of the CAGE window
						$flag++;
						print OUT join("\t",@a),"\n";
						print OUT1 join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},0),"\n";
					}
					elsif ($distL < 0) {
						my $min=abs($distL) < abs($distR) ? abs($distL) : abs($distR);
						if ($min <= $window) {
							$flag++;
							print OUT join("\t",@a),"\n";
							print OUT1 join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},(0-$min)),"\n";
						}
						else {
							$info=join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},(0-$min));
						}
					}
					elsif ($distL > 0) {
						my $min=$distL < $distR ? $distL : $distR;
						if ($min <= $window) {
							$flag++;
							print OUT join("\t",@a),"\n";
							print OUT1 join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},$min),"\n";
						}
						else {
							$info=join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},$min);
						}
					}
				}
			}
			$bin=$bin-1;
			if (($flag eq 0) and (exists $uniq{$a[2]}{$a[3]}{$bin})) {
				foreach my $pos(keys %{$uniq{$a[2]}{$a[3]}{$bin}}) {
					my @b=split("\t",$uniq{$a[2]}{$a[3]}{$bin}{$pos});
					my $distL=$a[4]-$b[0];
					my $distR=$a[4]-$b[1];
					if ($debug > 0) { print join("\t",$a[1],$a[2],$a[3],$a[4],$a[5],"bin",$bin,$distL,$distR),"\n"; }
					if ($distL * $distR <= 0) {
						# 5'end sit inside of the CAGE window
						$flag++;
						print OUT join("\t",@a),"\n";
						print OUT1 join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},0),"\n";
					}
					elsif ($distL < 0) {
						my $min=abs($distL) < abs($distR) ? abs($distL) : abs($distR);
						if ($min <= $window) {
							$flag++;
							print OUT join("\t",@a),"\n";
							print OUT1 join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},(0-$min)),"\n";
						}
						else {
							if ($info ne "") {
								my @tmp=split("\t",$info);
								if (abs($tmp[-1]) > $min) {
									$info=join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},(0-$min));
								}
							}
							else {
								$info=join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},(0-$min));
							}
						}
					}
					elsif ($distL > 0) {
						my $min=$distL < $distR ? $distL : $distR;
						if ($min <= $window) {
							$flag++;
							print OUT join("\t",@a),"\n";
							print OUT1 join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},$min),"\n";
						}
						else {
							$info=join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},$min);
							if ($info ne "") {
								my @tmp=split("\t",$info);
								if ($tmp[-1] > $min) {
									$info=join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},$min);
								}
							}
							else {
								$info=join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},$min);
							}
						}
					}
				}
			}
			$bin=$bin+2;
			if (($flag eq 0) and (exists $uniq{$a[2]}{$a[3]}{$bin})) {
				foreach my $pos(keys %{$uniq{$a[2]}{$a[3]}{$bin}}) {
					my @b=split("\t",$uniq{$a[2]}{$a[3]}{$bin}{$pos});
					my $distL=$a[4]-$b[0];
					my $distR=$a[4]-$b[1];
					if ($debug > 0) { print join("\t",$a[1],$a[2],$a[3],$a[4],$a[5],"bin",$bin,$distL,$distR),"\n"; }
					if ($distL * $distR <= 0) {
						# 5'end sit inside of the CAGE window
						$flag++;
						print OUT join("\t",@a),"\n";
						print OUT1 join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},0),"\n";
					}
					elsif ($distL < 0) {
						my $min=abs($distL) < abs($distR) ? abs($distL) : abs($distR);
						if ($min <= $window) {
							$flag++;
							print OUT join("\t",@a),"\n";
							print OUT1 join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},(0-$min)),"\n";
						}
						else {
							if ($info ne "") {
								my @tmp=split("\t",$info);
								if (abs($tmp[-1]) > $min) {
									$info=join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},(0-$min));
								}
							}
							else {
								$info=join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},(0-$min));
							}
						}
					}
					elsif ($distL > 0) {
						my $min=$distL < $distR ? $distL : $distR;
						if ($min <= $window) {
							$flag++;
							print OUT join("\t",@a),"\n";
							print OUT1 join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},$min),"\n";
						}
						else {
							$info=join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},$min);
							if ($info ne "") {
								my @tmp=split("\t",$info);
								if ($tmp[-1] > $min) {
									$info=join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},$min);
								}
							}
							else {
								$info=join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},$min);
							}
						}
					}
				}
			}
			if (($flag eq 0) and ($info ne "")) { print OUT1 $info,"\n"; }
		}
		else {
			my $bin=int($a[5]/$window);
			if ($debug > 0) { print "bin = ",$bin,"\n"; }
			my $flag=0;
			my $info="";
			if (exists $uniq{$a[2]}{$a[3]}{$bin}) {
				foreach my $pos(keys %{$uniq{$a[2]}{$a[3]}{$bin}}) {
					my @b=split("\t",$uniq{$a[2]}{$a[3]}{$bin}{$pos});
					my $distL=$a[5]-$b[0];
					my $distR=$a[5]-$b[1];
					if ($debug > 0) { print join("\t",$a[1],$a[2],$a[3],$a[4],$a[5],"bin",$bin,$distL,$distR),"\n"; }
					if ($distL * $distR <= 0) {
						# 5'end sit inside of the CAGE window
						$flag++;
						print OUT join("\t",@a),"\n";
						print OUT1 join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},0),"\n";
					}
					elsif ($distL < 0) {
						my $min=abs($distL) < abs($distR) ? abs($distL) : abs($distR);
						if ($min <= $window) {
							$flag++;
							print OUT join("\t",@a),"\n";
							print OUT1 join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},($min)),"\n";
						}
						else {
							$info=join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},($min));
						}
					}
					elsif ($distL > 0) {
						my $min=$distL < $distR ? $distL : $distR;
						if ($min <= $window) {
							$flag++;
							print OUT join("\t",@a),"\n";
							print OUT1 join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},(0-$min)),"\n";
						}
						else {
							$info=join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},(0-$min));
						}
					}
				}
			}
			$bin=$bin-1;
			if (($flag eq 0) and (exists $uniq{$a[2]}{$a[3]}{$bin})) {
				foreach my $pos(keys %{$uniq{$a[2]}{$a[3]}{$bin}}) {
					my @b=split("\t",$uniq{$a[2]}{$a[3]}{$bin}{$pos});
					my $distL=$a[5]-$b[0];
					my $distR=$a[5]-$b[1];
					if ($debug > 0) { print join("\t",$a[1],$a[2],$a[3],$a[4],$a[5],"bin",$bin,$distL,$distR),"\n"; }
					if ($distL * $distR <= 0) {
						# 5'end sit inside of the CAGE window
						$flag++;
						print OUT join("\t",@a),"\n";
						print OUT1 join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},0),"\n";
					}
					elsif ($distL < 0) {
						my $min=abs($distL) < abs($distR) ? abs($distL) : abs($distR);
						if ($min <= $window) {
							$flag++;
							print OUT join("\t",@a),"\n";
							print OUT1 join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},($min)),"\n";
						}
						else {
							if ($info ne "") {
								my @tmp=split("\t",$info);
								if (abs($tmp[-1]) > $min) {
									$info=join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},($min));
								}
							}
							else {
								$info=join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},($min));
							}
						}
					}
					elsif ($distL > 0) {
						my $min=$distL < $distR ? $distL : $distR;
						if ($min <= $window) {
							$flag++;
							print OUT join("\t",@a),"\n";
							print OUT1 join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},(0-$min)),"\n";
						}
						else {
							$info=join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},(0-$min));
							if ($info ne "") {
								my @tmp=split("\t",$info);
								if ($tmp[-1] > $min) {
									$info=join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},(0-$min));
								}
							}
							else {
								$info=join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},(0-$min));
							}
						}
					}
				}
			}
			$bin=$bin+2;
			if (($flag eq 0) and (exists $uniq{$a[2]}{$a[3]}{$bin})) {
				foreach my $pos(keys %{$uniq{$a[2]}{$a[3]}{$bin}}) {
					my @b=split("\t",$uniq{$a[2]}{$a[3]}{$bin}{$pos});
					my $distL=$a[5]-$b[0];
					my $distR=$a[5]-$b[1];
					if ($debug > 0) { print join("\t",$a[1],$a[2],$a[3],$a[4],$a[5],"bin",$bin,$distL,$distR),"\n"; }
					if ($distL * $distR <= 0) {
						# 5'end sit inside of the CAGE window
						$flag++;
						print OUT join("\t",@a),"\n";
						print OUT1 join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},0),"\n";
					}
					elsif ($distL < 0) {
						my $min=abs($distL) < abs($distR) ? abs($distL) : abs($distR);
						if ($min <= $window) {
							$flag++;
							print OUT join("\t",@a),"\n";
							print OUT1 join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},($min)),"\n";
						}
						else {
							if ($info ne "") {
								my @tmp=split("\t",$info);
								if (abs($tmp[-1]) > $min) {
									$info=join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},($min));
								}
							}
							else {
								$info=join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},($min));
							}
						}
					}
					elsif ($distL > 0) {
						my $min=$distL < $distR ? $distL : $distR;
						if ($min <= $window) {
							$flag++;
							print OUT join("\t",@a),"\n";
							print OUT1 join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},(0-$min)),"\n";
						}
						else {
							$info=join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},(0-$min));
							if ($info ne "") {
								my @tmp=split("\t",$info);
								if ($tmp[-1] > $min) {
									$info=join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},(0-$min));
								}
							}
							else {
								$info=join("\t",$a[0],$a[1],$a[2],$a[3],$uniq{$a[2]}{$a[3]}{$bin}{$pos},(0-$min));
							}
						}
					}
				}
			}
			if (($flag eq 0) and ($info ne "")) { print OUT1 $info,"\n"; }
		}
	}
	else {
		
	}
}
