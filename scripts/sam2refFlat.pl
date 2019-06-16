#!/usr/bin/perl -w
use strict;
# 2014-05-12 by Arthur
die "Usage: $0   \"input_sam\"   \"output_refFlat\"    \(optional\)\"Chromosome_name_replacement\"   \(DEBUG==1\)" if (@ARGV < 2);
my $filein1=$ARGV[0];
my $fileout=$ARGV[1];
my $chromo="";
if (scalar(@ARGV) > 2) {$chromo=$ARGV[2];}	# /home/arthur/project/Erin/PB_rat/0_raw_gmap/idlist
my $debug=0;
if (scalar(@ARGV) > 3) {$debug=$ARGV[3];}
my %Chromo;
if ($chromo ne "") {
	open IN, $chromo;
	while(<IN>) {
		chomp;
		my @a=split("\t",$_);	# gi|380690196|gb|CM000072.4|     chr1
		$Chromo{$a[0]}=$a[1];
	}
	close IN;
}

sub dec2bin {
    my $str = unpack("B32", pack("N", shift));
    #$str =~ s/^0+(?=\d)//; # otherwise you'll get leading zeros return $str;
    my $decode = substr($str,-11);
    return $decode;
}

sub bin2dec {
    return unpack("N", pack("B32", substr("0" x 32 . shift, -32)));
}

open(IN1, $filein1) or die "Cannot open input_sam file : $filein1";
open(OUT, ">".$fileout) or die "Cannot open output file : $fileout";
while (<IN1>) {
    chomp;
    if ((m/^@/) or (m/^#/)) { next; }
    my @a=split("\t",$_);
    # 0                                                     1    2                              3            4      5                                                    6    7    8    9      10   11      
    # RNN_LD2.104_ccs_5502_2048_0_2048_2048_SUB0_pTNp_17    0    gi|380690193|gb|CM000075.4|    192723046    255    6S447M2042N121M606631N141M109471N125M5559N1129M2S    *    0    0    seq    *    NH:i:1  HI:i:1  AS:i:1963    nM:i:0
	next if ($a[1] eq 4);
    my $cFlag=dec2bin($a[1]);
    my @aFlag=split('',$cFlag);
    my $strand="+";
    if ($aFlag[6] eq 1) {$strand="-";}
	my $chr=$a[2];
	if (($chromo ne "") and (exists $Chromo{$a[2]})) { $chr=$Chromo{$a[2]}; }
	my @Start;
	my @End;
	my @CIGAR_op=($a[5]=~m/[MSNID]/g); 
    my @CIGAR_va=($a[5]=~m/\d+/g);
    my $start=$a[3]-1;
	my $exon=0;
    my $Nr=scalar(@CIGAR_op);
	for(my $i=0; $i<$Nr; $i++){
		# leading and ending mark "S" are ignored; 
	    if ($CIGAR_op[$i] eq "M") {
			$Start[$exon]=$start;
			$start+=$CIGAR_va[$i];
			$End[$exon]=$start;
		}
	    elsif ($CIGAR_op[$i] eq "D") {
			$start+=$CIGAR_va[$i];
			$End[$exon]=$start;
		}
	    elsif ($CIGAR_op[$i] eq "N") {
			$exon++;
			$start+=$CIGAR_va[$i];
	    }
	}
	$exon++;
	print OUT join("\t",$a[0],$a[0],$chr,$strand,$Start[0],$End[-1],$Start[0],$Start[0],$exon,join(",",@Start),join(",",@End)),"\n";
}


