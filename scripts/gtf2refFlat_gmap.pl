#!/usr/bin/perl -w
use strict;
# 2014-02-16 by Arthur
die "Usage: $0   \"input_gtf\"    \"input_len\"    \"output_refFlat\"   \(optional\)\"use_best==1, set to 0 to discard multiple alignment\"  \(optional\)\"pacbio==1\"   \(optional\)\"Chromosome_name_replacement\"   \(DEBUG==1\)" if (@ARGV < 3);
my $filein=$ARGV[0];
my $filelen=$ARGV[1];
my $fileout=$ARGV[2];
my $useBest=1;
if (scalar(@ARGV) > 3) {$useBest=$ARGV[3];}
my $pacbio=0;
if (scalar(@ARGV) > 4) {$pacbio=$ARGV[4];}
my $chromo="";
if (scalar(@ARGV) > 5) {$chromo=$ARGV[5];}	# /home/arthur/project/Erin/PB_rat/0_raw_gmap/idlist
my $debug=0;
if (scalar(@ARGV) > 6) {$debug=$ARGV[6];}
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

my %LEN;
open(IN, $filelen) or die "cannot open length_file : $filelen";
while (<IN>) {
	chomp;
	next if (m/^#/);
	s/^>//;
	my @a=split("\t",$_);
	$LEN{$a[0]}=$a[1];
}
close IN;


open IN,$filein;
open OUT,">".$fileout;
open OUT1,">".$fileout.".cov";
open MUL,">".$fileout.".mul";
open ERR,">".$fileout.".err";
my %Gene;
my %COV;	# coverage
my %IDEN;	# identity
my %Covered;# percentage of nts that is covered by sequence
my %Matches;# number of matched nts
my %Multi;	# multiple hits (if and only if they are NOT on the same_chromosome + same_strand)
my %Exon;
my %CDS;
my %Read;	# align-able border of the read
my %WELL;	# pacbio_sequencing_well				# $WELL{$Tid}=$WellId;
my %BadW;	# if one subread from a well is problematic, mark all subreads from that very well	# $BadW{$WellId}=1;
my %exonQua;# the alignment quality of each exons, integer ranging in [0,100]
my $newmrna=0;
while(<IN>) {
    next if(/^#/);
    chomp;
    my $info=$_;
    my @a=split("\t",$info);
	if (($a[0]=~m/^chr/) or ($a[0]=~m/^Chr/)){}
	elsif (exists $Chromo{$a[0]}) {my $t=$Chromo{$a[0]}; $a[0]=$t;}
	$a[0]=~s/\_/\|/g;	# this is to "rescue" chr named as "chrY_KL568141v1_random"
	# else {next;}
	my $Tid="na";	# Transcript
	my $Gid="na";	# Gene
    if ($info=~m/ID=([\w\d\.\|]+).mrna/) {$Tid=$1;}
	else {next;}
	if (($pacbio eq 1) and ($a[2] eq "mRNA")) {
		my @tmp=split(/\_/,$Tid);
		my $flag=0;
		for(my $i=0; $i<scalar(@tmp); $i++) {
			if (($tmp[$i] eq "ccs") or ($tmp[$i] eq "sub")) {
				$flag=$i;
				last;
			}
		}
		my $WellId=$tmp[0];
		for(my $i=1; $i<=($flag+2); $i++) {$WellId=$WellId."_".$tmp[$i];}
		$WELL{$Tid}=$WellId;
	}
    if ($info=~m/ID=(Isogroup_\d+)_\d+_\d+_\d+/) { $Gid=$1; }
    elsif (/gene_name\s\"(\S+)\";/) {$Gid=$1;}
    else {$Gid=$Tid;}
	if ($debug eq 1) { print $Tid,"\t",$Gid,"\n"; }
    $a[3]-=1;
	if ($a[2] eq "mRNA") {
		$newmrna=0;
		my @b=split(/\;/,$a[8]);
		# always keep the one with most matches, and reset the $Read{}, $CDS{}, $Exon{} if exists
		my $tmp_mat=0;
		my $tmp_cov=0;
		my $tmp_ide=0;
		for(@b) {
			my @c=split(/\=/,$_);
			if ($c[0] eq "coverage") {$tmp_cov=$c[1];} 
			if ($c[0] eq "identity") {$tmp_ide=$c[1];}
			if ($c[0] eq "matches") {$tmp_mat=$c[1];}
		}
		if (exists $Matches{$Tid}) {
            if ($Matches{$Tid} < $tmp_mat) {
                $COV{$Tid}=$tmp_cov;
				$IDEN{$Tid}=$tmp_ide;
				$Matches{$Tid}=$tmp_mat;
				delete $CDS{$Tid};
				delete $Exon{$Tid};
				delete $Read{$Tid};
				delete $Covered{$Tid};
				delete $Gene{$Tid};
				delete $exonQua{$Tid};
				if ($useBest eq 0) {
                    $Multi{$Tid}=1;
					if ($pacbio eq 1) { $BadW{$WELL{$Tid}}=1; }
					print MUL join("\t",$Tid,$Gid,$Gene{$Tid},$Gid."\_\_\_\_".$a[0]."\_\_\_\_".$a[6]."\_\_\_\_".$b[3]),"\n";
                }
            }
			else { $newmrna=1; }
        }
		else {
			$COV{$Tid}=$tmp_cov;
			$IDEN{$Tid}=$tmp_ide;
			$Matches{$Tid}=$tmp_mat;
		}
	}
	elsif (($a[2] eq "exon") and ($newmrna eq 0)) {
		if (!exists $Gene{$Tid}) {
			#							  chr			   strand_on_genome strand_on_sequence
			my @b=split(" ",$a[8]);
			$Gene{$Tid} = $Gid."\_\_\_\_".$a[0]."\_\_\_\_".$a[6]."\_\_\_\_".$b[3];
		}
		my $str=join("_",$a[0],$a[6],$a[3],$a[4]);
		if (exists $Exon{$Tid} ) { $Exon{$Tid}=$Exon{$Tid}."\t".$str; }
		else { $Exon{$Tid}=$str; }
		my @b=split(" ",$a[8]);
		my $str2=join("\t",$b[1],$b[2],$b[3]);
		if (exists $Read{$Tid}) {
			my @c=split("\t",$Read{$Tid});
			my @d;
			push @d, $c[0];
			push @d, $c[1];
			push @d, $b[1];
			push @d, $b[2];
			my @sd=sort{$a <=> $b} @d;
			$Read{$Tid}=join("\t",$sd[0],$sd[3],$c[2]);
			$exonQua{$Tid}=$exonQua{$Tid}."\t".$a[5];
		}
		else {
			$Read{$Tid}=$str2;
			$exonQua{$Tid}=$a[5];
		}
		$Covered{$Tid}+=$a[4]-$a[3];
		
	}
	elsif (($a[2] eq "CDS") and ($newmrna eq 0)) {
	    my $str=join("_",$a[0],$a[6],$a[3],$a[4]);
	    if (exists $CDS{$Tid} ) { $CDS{$Tid} = $CDS{$Tid}."\t".$str; }
		else { $CDS{$Tid} = $str; }
	}
	
}
# close MUL;

foreach my $Tid ( sort {$Gene{$a} cmp $Gene{$b} } keys %Gene ) {
	if (exists $Multi{$Tid}) {next;}
	if (($pacbio eq 1) and (exists $BadW{$WELL{$Tid}})) {
		print MUL $Tid,"\n";
		next;
	}
	my @a=split("\t",$Read{$Tid});
	my $cov=sprintf("%.1f",100*($a[1] - $a[0])/$LEN{$Tid});
	print OUT1 join("\t",$Tid,$cov,$IDEN{$Tid},$Read{$Tid},$LEN{$Tid},$Covered{$Tid},$COV{$Tid},$exonQua{$Tid}),"\n";
	my ($cdsStart, $cdsEnd, $exonStarts, $exonEnds);
	my $transID=$Tid;
	@a=split("\_\_\_\_",$Gene{$Tid});
	my $chr=$a[1];
	my $strand=$a[2];
	# process exons
	my @starts;
	my @ends;
	my @b=split("\t",$Exon{$Tid});
	my $Nrb=scalar(@b);
	for(my $i=0; $i<$Nrb; $i++) {
		my @c=split(/\_/,$b[$i]);
		if ($chr ne $c[0]) {print ERR join("\t","Exon_Chr",$Tid,$chr,$c[0]),"\n";}
		if ($strand ne $c[1]) {print ERR join("\t","Exon_strand",$Tid,$strand,$c[1]),"\n";}
		push @starts, $c[2];
		push @ends, $c[3];
	}
	my @Sstart=sort {$a <=> $b} @starts;
	my @Send=sort {$a <=> $b} @ends;
	my $start=$Sstart[0];
	my $end=$Send[$Nrb - 1];
	my $exon_starts="";
	my $exon_ends="";
	for(my $i=0; $i<$Nrb; $i++) {
		$exon_starts=$exon_starts.$Sstart[$i].",";
		$exon_ends=$exon_ends.$Send[$i].",";
	}
	# process CDS
	my $cds_start=$start;
	my $cds_end=$start;
	if (exists $CDS{$Tid}) {
		my @starts;
		my @ends;
		my @d=split("\t",$CDS{$Tid});
		my $Nrd=scalar(@d);
		for(my $i=0; $i<$Nrd; $i++) {
			my @e=split(/\_/,$d[$i]);
			if ($chr ne $e[0]) {print ERR join("\t","CDS_Chr",$Tid,$chr,$e[0]),"\n";}
			if ($strand ne $e[1]) {print ERR join("\t","CDS_Strand",$Tid,$chr,$e[1]),"\n";}
			push @starts, $e[2];
			push @ends, $e[3];
		}
		my @CDS_start=sort {$a <=> $b} @starts;
		my @CDS_end=sort {$a <=> $b} @ends;
		$cds_start=$CDS_start[0];
		$cds_end=$CDS_end[$Nrd - 1];
	}
	# output
	if ($Nrb > 0) {
		$a[0]=~s/\|/\_/g;
		print OUT join("\t",$a[0],$Tid,$chr,$strand,$start,$end,$cds_start,$cds_end,$Nrb,$exon_starts,$exon_ends),"\n";
	}
}
close IN;
close OUT;
close MUL;
my $fileu1=$fileout.".mul";
my $fileu2=$fileout.".multi";
my $command="uniq $fileu1 > $fileu2";
system($command);
$command="rm -f $fileu1";
system($command);
