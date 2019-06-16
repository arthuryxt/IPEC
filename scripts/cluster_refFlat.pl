#!/usr/bin/perl -w
use strict;
# 2014-02-16 by Arthur
die "Usage: $0   \"input_refFlat\"    \"output_refFlat.basename\"   \"seqlenth\"   \"strand_info\(default=\"no\"\)\"   \"coverage_file\(default=\"no\"\)\"    \"annotated_junction_bed6_bed12\"\(default=\"no\"\)   \"Illumina_junction_bed6_bed12\"\(default=\"no\"\)    \"CAGE_bed6\"\(default=\"no\"\)    \"min_junction_support==2\"   \"min_intron_len==50\"  \"max_intron_len==100000\"   \"min_isof_len==300\"  \(force_illumina_support==0\)   \(DEBUG==1\)" if (@ARGV < 2);
my $filein1=$ARGV[0];
my $fileout=$ARGV[1];
my $seqlength=$ARGV[2];
my $strandinfo="no";
if (scalar(@ARGV) > 3) { $strandinfo=$ARGV[3]; }
my $coverage="no";
if (scalar(@ARGV) > 4) { $coverage=$ARGV[4]; }
my $support="no";	# any support of junction
if (scalar(@ARGV) > 5) { $support=$ARGV[5]; }
my $knownjunc="no";	# Illumina_support for screening
if (scalar(@ARGV) > 6) { $knownjunc=$ARGV[6]; }
my $knowncage="no";
if (scalar(@ARGV) > 7) { $knowncage=$ARGV[7]; }
my $min_junc=2;
if (scalar(@ARGV) > 8) { $min_junc=$ARGV[8]; }
my $min_intron=50;
if (scalar(@ARGV) > 9) { $min_intron=$ARGV[9]; }
my $max_intron_len=100000;	# longest_gene_hg19 = 4349516; longest_intron_hg19 = 91666; longest_intron_mm9 = 52269;
if (scalar(@ARGV) > 10) { $max_intron_len=$ARGV[10]; }
my $min_isof_len=300;
if (scalar(@ARGV) > 11) { $min_isof_len=$ARGV[11]; }
my $force_illumina=0;
if (scalar(@ARGV) > 12) { $force_illumina=$ARGV[12]; }
my $debug=998;
if (scalar(@ARGV) > 13) { $debug=$ARGV[13]; open(OUTdebug,">".$fileout.".debug");  }

my $min_3UTR_diff=100;
my $fuzzy=10;				# fuzziness of gmap exon-border
my $CovCutOff=85;			# at least 85% of the sequence can be covered by RNA-Seq reads
my $aln_min_cov=0.9;		# at least 90% of the sequence can be aligned to the genome

print join("\t","input_refFlat :",$filein1),"\n";
print join("\t","output file   :",$fileout),"\n";
print join("\t","seqlength     :",$seqlength),"\n";
print join("\t","strandinfo    :",$strandinfo),"\n";
print join("\t","coverage      :",$coverage),"\n";
print join("\t","anno_junc     :",$support),"\n";
print join("\t","RNA-Seq_junc  :",$knownjunc),"\n";
print join("\t","min_junc_len  :",$min_junc),"\n";
print join("\t","min_intron_len:",$min_intron),"\n";
print join("\t","max_intron_len:",$max_intron_len),"\n";
print join("\t","min_isof_len  :",$min_isof_len),"\n";
print join("\t","force_illumina:",$force_illumina),"\n";
print join("\t","debug         :",$debug),"\n";

# store the coverage information of each read
my %Coverage;
if ($coverage ne "no") {
	open(IN, $coverage) or die "Cannot open coverage file $coverage ! Type no to ingore coverage_info.";
	while (<IN>) {
		chomp;
		my @a=split("\t",$_);
		$Coverage{$a[0]}=join("\t",@a);
	}
	close IN;
}

# store the information of each read : strandness and category
my %Strand;
if ($strandinfo ne "no") {
	open(IN, $strandinfo) or die "Cannot open strand_info file $strandinfo ! Type no to ingore input_info.";
	while (<IN>) {
		chomp;
		my @a=split("\t",$_);
		# Readid	strand	Category : use strand_info ONLY for single-exon read;	use Class0 ONLY for determine TSS
		$Strand{$a[0]}=join("\t",@a);
	}
	close IN;
}

# store the support_junctions from Annotations
my %SupportJunc;
my %SupportExonBorder;
if ($support ne "no") {
	open(IN, $support) or die "Cannot open bed6- or bed12-format known_junction file $support ! Type no to ignore known_junction.";
	while (<IN>) {
		chomp;
		if (m/track/) {next;}
		if (m/^Junction/) {next;}
		my @a=split("\t",$_);
		# 1		312360	315196	ENSRNOT00000044669	1	+
		# 1		4482658 4483189 JUNC00000001		1	-       4482658 4483189 255,0,0 2       91,9    0,522
		if ($a[0]!~m/^chr/) {my $chr="chr".$a[0]; $a[0]=$chr; }
		my $left=$a[1];
		my $right=$a[2];
		if (scalar(@a) eq 12) {
			my @b=split(/\,/,$a[10]);
			$left=$a[1]+$b[0];
			$right=$a[2]-$b[1];
		}
		if (($right - $left) < $min_intron) { next; }
		if (($right - $left) > $max_intron_len) { next; }
		$SupportJunc{$a[0]}{$a[5]}{$left."\t".$right}=1;
		$SupportExonBorder{$a[0]}{$a[5]}{$left}=1;
		$SupportExonBorder{$a[0]}{$a[5]}{$right}=1;
	}
	close IN;
}

# store the known_junctions from Illumina
my %KnownJunc;
my %KnownExonBorder;
if ($knownjunc ne "no") {
	open(IN, $knownjunc) or die "Cannot open bed6- or bed12-format known_junction file $knownjunc ! Type no to ignore known_junction.";
	while (<IN>) {
		chomp;
		if (m/track/) {next;}
		my @a=split("\t",$_);
		# 1       4482658 4483189 JUNC00000001    1       -       4482658 4483189 255,0,0 2       91,9    0,522
		if ($a[0]!~m/^chr/) {my $chr="chr".$a[0]; $a[0]=$chr; }
		my @b=split(/\,/,$a[10]);
		my $left=$a[1];
		my $right=$a[2];
		if (scalar(@a) eq 12) {
			$left=$a[1]+$b[0];
			$right=$a[2]-$b[1];
		}
		if (($right - $left) < $min_intron) { next; }
		my $bin=int($right/1000);
		$KnownJunc{$a[0]}{$a[5]}{$bin}{$left."\t".$right}+=$a[4];
		$KnownExonBorder{$a[0]}{$a[5]}{$left}+=$a[4];
		$KnownExonBorder{$a[0]}{$a[5]}{$right}+=$a[4];
	}
	close IN;
}
if ($debug > 0){
    open(OUTKJ, ">KnownJunc");
	foreach my $chr (sort keys %KnownExonBorder) {
		foreach my $strand (sort keys %{$KnownExonBorder{$chr}}) {
			foreach my $pos (sort keys %{$KnownExonBorder{$chr}{$strand}}) {
				print OUTKJ join("\t",$chr,$strand,$pos,$KnownExonBorder{$chr}{$strand}{$pos}),"\n";
			}
		}
	}
	close OUTKJ;
}


# store the known_Cage
my %KnownCage;
if ($knowncage ne "no") {
	open(IN, $knowncage) or die "Cannot open bed6-format CAGE file $knowncage ! Type no to ignore Cage.";
	while (<IN>) {
		chomp;
		my @a=split("\t",$_);
		# chr1	start	end	name	score	strand
		my $bin=int($a[1]/1000);
		$KnownCage{$a[0]}{$a[5]}{$bin}{$a[1]}=join("\t",@a);
	}
	close IN;
}

# store the length of raw sequence
my %rawseqlen;
open(IN,$seqlength) or die "Cannot open seqlength file $seqlength";
while(<IN>) {
	chomp;
	my @a=split("\t",$_);
	$rawseqlen{$a[0]}=$a[1];
}
close IN;

# read refFlat one by one, record by chr_strand_start, make gene_loci chr_strand_start_end, store exact_junc(fix fuzziness later.)
my %Gene_loci;
my %Junc;		# all junctions in a exact way
my %Read;		# info for every read
my %Label;		# gene_loci for every reads
my %Len;
my %Single;		# single-exon-read
my %RareRead;
my %Truncate;
open(IN, $filein1) or die "Cannot open input_refFlat $filein1";
open(OUT0,">".$fileout.".loci");
open(OUT1,">".$fileout.".cluster");	
open(OUT2,">".$fileout.".short");
open(OUT3,">".$fileout.".long");
open(OUT4,">".$fileout.".truncated");
open(OUT5,">".$fileout.".rare");
open(OUT6,">".$fileout.".single");
open(OUT7,">".$fileout.".lowcov");
open(OUT8,">".$fileout.".Junc");
open(OUT9,">".$fileout.".mExon");
open(OUT10, ">".$fileout.".negativeExonLen");
while (<IN>) {
	chomp;
	my @a=split("\t",$_);
	if ($a[2]!~m/^chr/) {my $chr="chr".$a[2]; $a[2]=$chr; }
	## discard read if of low coverage
	if ($coverage ne "no") {
		if (exists $Coverage{$a[0]}) {
			my @b=split("\t",$Coverage{$a[0]});
			if ($b[1] < $CovCutOff) {
				print OUT7 join("\t",@a),"\t",join("\t",@b),"\n";
				next;
			}
		}
		else {
			next;
		}
	}
	## single-exon-read-without-known-strandness
	if ($a[8] eq 1) {
		if (($strandinfo ne "no") and (exists $Strand{$a[0]})) {
			my @b=split("\t",$Strand{$a[0]});
			$a[3]=$b[1];
		}
		else { $Single{$a[0]}=1; }
	}
	# check strand info # this would lead to chaos due to pseudogene
	#if (($strandinfo ne "no") and (exists $Strand{$a[0]})) {
	#	my @b=split("\t",$Strand{$a[0]});
	#	$a[3]=$b[1];
	#}
	if($a[3] eq "na") {next;}
	if($a[3] eq "+") {$a[6] = $a[4]; $a[7] = $a[4];}
	if($a[3] eq "-") {$a[6] = $a[5]; $a[7] = $a[5];}
	my @Start=split(/\,/,$a[9]);
	my @End=split(/\,/,$a[10]);
	my $len=0;
	for(my $i=0; $i<$a[8]; $i++) { $len = $len + $End[$i] - $Start[$i]; }
	$Len{$a[0]}=$len;
	my $tmpalnpct=sprintf("%.3f",($len/$rawseqlen{$a[0]}));
	if($len < $min_isof_len) {print OUT2 join("\t",@a),"\t",$len,"\t",$rawseqlen{$a[0]},"\t",$tmpalnpct,"\n"; next; }
	if($tmpalnpct < $aln_min_cov) {print OUT2 join("\t",@a),"\t",$len,"\t",$rawseqlen{$a[0]},"\t",$tmpalnpct,"\n"; next; }
	if($tmpalnpct > 1.05) {print OUT2 join("\t",@a),"\t",$len,"\t",$rawseqlen{$a[0]},"\t",$tmpalnpct,"\n"; next; }
	my $MIntron=-1;
	for(my $i=1; $i<$a[8]; $i++) { $len=abs($Start[$i] - $End[$i-1]); if ($len > $MIntron){$MIntron=$len};}
	if ($MIntron > $max_intron_len) { print OUT3 join("\t",@a),"\t",$MIntron,"\n"; next; }
	## store read_info
	$Read{$a[0]}=join("\t",@a);
	if (($a[8] eq 1) and (exists $Single{$a[0]})) { $a[3]="na"; }
	## check gene_loci
	# my $start=100*(int($a[4]/100));
	# my $end=100*(int($a[5]/100))+100;
	# if (exists $Gene_loci{$a[2]}) {
	# 	if (exists $Gene_loci{$a[2]}{$a[3]}) {
	# 		my $flag=0;
	# 		my %conflict;
	# 		foreach my $pos(sort{$a <=> $b} keys %{$Gene_loci{$a[2]}{$a[3]}}) {
	# 			if ($Gene_loci{$a[2]}{$a[3]}{$pos} < $start) {	}
	# 			elsif ($end < $pos) { last; }
	# 			else {
	# 				$conflict{$pos}=$Gene_loci{$a[2]}{$a[3]}{$pos};
	# 				$flag++;
	# 			}
	# 		}
	# 		if ($flag eq 0) { $Gene_loci{$a[2]}{$a[3]}{$start}=$end; }
	# 		else {
	# 			my @Pos;
	# 			push @Pos,$start;
	# 			push @Pos,$end;
	# 			foreach my $pos (keys %conflict) {
	# 				push @Pos,$pos;
	# 				push @Pos,$conflict{$pos};
	# 				delete $Gene_loci{$a[2]}{$a[3]}{$pos};
	# 			}
	# 			my @POS=sort{$a <=> $b} @Pos;
	# 			$Gene_loci{$a[2]}{$a[3]}{$POS[0]}=$POS[2*$flag+1];
	# 		}
	# 	}
	# 	else { $Gene_loci{$a[2]}{$a[3]}{$start}=$end; }
	# }
	# else { $Gene_loci{$a[2]}{$a[3]}{$start}=$end; }
	# if ($debug <= 50) {print OUTdebug join("\t",$a[2],$a[3],$start,$end,$Read{$a[0]}),"\n";}
	## store fuzzy_junc
	for(my $i=1; $i<$a[8]; $i++) {
		# 1       4482658 4483189 JUNC00000001    1       -       4482658 4483189 255,0,0 2       91,9    0,522
		if ($debug <= 10) { print OUTdebug join("\t",$a[2],$End[$i-1],$Start[$i],"junc",$a[3],$End[$i-1],$Start[$i],join(",",255,0,0),2,join(",",0,0),join(",",0,($Start[$i] - $End[$i-1]))),"\n"; }
		$Junc{$a[2]."_".$a[3]."_".$End[$i-1]."_".$Start[$i]}++;
	}
}
close IN;
close OUT2;

my %NewJunc;
my %ModJunc;
foreach my $id (sort{$Junc{$b} <=> $Junc{$a}} keys %Junc) {
	my @a=split(/\_/,$id);
	#my $bin=int($a[2]/1000);
	# check if is Known
	if ((exists $KnownExonBorder{$a[0]}{$a[1]}{$a[2]}) and ($KnownExonBorder{$a[0]}{$a[1]}{$a[3]})) {
		$NewJunc{$a[0]}{$a[1]}{$a[2]."\t".$a[3]}=$Junc{$id};
	}
	# check if near Known
	else {
		my $flag="";
		for(my $i=0-$fuzzy; $i<=$fuzzy; $i++) {
			for(my $j=0-$fuzzy; $j<=$fuzzy; $j++) {
				if (exists $KnownJunc{$a[0]}{$a[1]}{($a[2]+$i)."\t".($a[3]+$j)}) {
					$flag=$i."\t".$j;
					last;
				}
			}
			if ($flag ne "") { last; }
		}
		if ($flag eq "") {
			# check if near recorded
			for(my $i=0-$fuzzy; $i<=$fuzzy; $i++) {
				for(my $j=0-$fuzzy; $j<=$fuzzy; $j++) {
					if (exists $NewJunc{$a[0]}{$a[1]}{($a[2]+$i)."\t".($a[3]+$j)}) {
						$flag=$i."\t".$j;
						last;
					}
				}
				if ($flag ne "") { last; }
			}
			if ($flag eq "") {
				$NewJunc{$a[0]}{$a[1]}{$a[2]."\t".$a[3]}=$Junc{$id};
			}
			else {
				my @t=split("\t",$flag);
				$ModJunc{$id}=join("_",$a[0],$a[1],($a[2]+$t[0]),($a[3]+$t[1]));
				$NewJunc{$a[0]}{$a[1]}{($a[2]+$t[0])."\t".($a[3]+$t[1])}+=$Junc{$id};
			}
		}
		else {
			my @t=split("\t",$flag);
			$ModJunc{$id}=join("_",$a[0],$a[1],($a[2]+$t[0]),($a[3]+$t[1]));
			$NewJunc{$a[0]}{$a[1]}{($a[2]+$t[0])."\t".($a[3]+$t[1])}+=$Junc{$id};
		}
	}	
}

# report Junction expression
my %JuncExp;
my %Jsample;
my %Exon;
foreach my $id (keys %Read) {
	if (exists $Single{$id}) { next; }
	else {
		my @aid=split(/\_/,$id);
		$Jsample{$aid[0]}=1;
		my @a=split("\t",$Read{$id});
		if ($a[8] < 2) {next;}
		my @Start=split(/\,/,$a[9]);
		my @End=split(/\,/,$a[10]);
		for(my $k=1; $k<$a[8]; $k++) {
			my $Jid=$a[2]."_".$a[3]."_".$End[$k-1]."_".$Start[$k];
			my $Jstart=$End[$k-1];
			my $Jend=$Start[$k];
			if (exists $ModJunc{$Jid}) {
				my @b=split(/\_/,$ModJunc{$Jid});
				$Jstart=$b[-2];
				$Jend=$b[-1];
			}
			$Jid=$a[2]."_".$a[3]."_".$Jstart."_".$Jend;
			$JuncExp{$Jid}{$aid[0]}++;
		}
		# first exon, check only the End
		my $k=0;
		my $Jid=$a[2]."_".$a[3]."_".$End[$k]."_".$Start[$k+1];
		my $Jstart=$End[$k];
		my $Jend=$Start[$k+1];
		if (exists $ModJunc{$Jid}) {
			my @b=split(/\_/,$ModJunc{$Jid});
			$Jstart=$b[-2];
			$Jend=$b[-1];
		}
		my $Eid=$a[2]."_".$a[3]."_".$Start[$k]."_".$Jstart."_tss";
		if ($a[3] eq "-") { $Eid=$a[2]."_".$a[3]."_".$Start[$k]."_".$Jstart."_end"; }
		$Exon{$Eid}{$aid[0]}++;
		# middle exons, check both Start and End
		for(my $k=1; $k<($a[8]-1); $k++) {
			my $Jid_1=$a[2]."_".$a[3]."_".$End[$k-1]."_".$Start[$k];
			my $Jstart_1=$End[$k-1];
			my $Jend_1=$Start[$k];
			if (exists $ModJunc{$Jid_1}) {
				my @b=split(/\_/,$ModJunc{$Jid_1});
				$Jstart_1=$b[-2];
				$Jend_1=$b[-1];
			}
			my $Jid_2=$a[2]."_".$a[3]."_".$End[$k]."_".$Start[$k+1];
			my $Jstart_2=$End[$k];
			my $Jend_2=$Start[$k+1];
			if (exists $ModJunc{$Jid_2}) {
				my @b=split(/\_/,$ModJunc{$Jid_2});
				$Jstart_2=$b[-2];
				$Jend_2=$b[-1];
			}
			my $Eid=$a[2]."_".$a[3]."_".$Jend_1."_".$Jstart_2."_mid";
			$Exon{$Eid}{$aid[0]}++;
		}
		# last exon,check only the Start
		$k=$a[8]-1;
		$Jid=$a[2]."_".$a[3]."_".$End[$k-1]."_".$Start[$k];
		$Jstart=$End[$k-1];
		$Jend=$Start[$k];
		if (exists $ModJunc{$Jid}) {
			my @b=split(/\_/,$ModJunc{$Jid});
			$Jstart=$b[-2];
			$Jend=$b[-1];
		}
		$Eid=$a[2]."_".$a[3]."_".$Jend."_".$End[$k]."_end";
		if ($a[3] eq "-") { $Eid=$a[2]."_".$a[3]."_".$Start[$k]."_".$Jstart."_tss"; }
		$Exon{$Eid}{$aid[0]}++;
	}
}
my $jreport="Junction\tLeft\tRight\tSum";
foreach my $sample (sort keys %Jsample) { $jreport=$jreport."\t".$sample; }
print OUT8 $jreport,"\n";
foreach my $Jid (sort keys %JuncExp) {
	my @Ja=split(/\_/,$Jid);
	$jreport=$Jid;
	if (exists $KnownExonBorder{$Ja[0]}{$Ja[1]}{$Ja[2]}) { $jreport=$jreport."\t".$KnownExonBorder{$Ja[0]}{$Ja[1]}{$Ja[2]}; }
	else { $jreport=$jreport."\t0"; }
	if (exists $KnownExonBorder{$Ja[0]}{$Ja[1]}{$Ja[3]}) { $jreport=$jreport."\t".$KnownExonBorder{$Ja[0]}{$Ja[1]}{$Ja[3]}; }
	else { $jreport=$jreport."\t0"; }
	$jreport=$jreport."\t".$NewJunc{$Ja[0]}{$Ja[1]}{$Ja[2]."\t".$Ja[3]};
	foreach my $sample (sort keys %Jsample) {
		if (exists $JuncExp{$Jid}{$sample}) { $jreport=$jreport."\t".$JuncExp{$Jid}{$sample} }
		else { $jreport=$jreport."\t0"; }
	}
	print OUT8 $jreport,"\n";
}
close OUT8;
$jreport="Exon";
foreach my $sample (sort keys %Jsample) { $jreport=$jreport."\t".$sample; }
print OUT9 $jreport,"\n";
foreach my $Jid (sort keys %Exon) {
	my @Ja=split(/\_/,$Jid);
	$jreport=$Jid;
	foreach my $sample (sort keys %Jsample) {
		if (exists $Exon{$Jid}{$sample}) { $jreport=$jreport."\t".$Exon{$Jid}{$sample} }
		else { $jreport=$jreport."\t0"; }
	}
	print OUT9 $jreport,"\n";
}
close OUT9;

# make Gene_loci with multi-exon reads with good junctions
foreach my $id (keys %Read) {
	if (exists $Single{$id}) { next; }
	else {
		my @a=split("\t",$Read{$id});
		if ($a[8] < 2) {next;}
		my @Start=split(/\,/,$a[9]);
		my @End=split(/\,/,$a[10]);
		# check if there's no poor-junctions
		my $RARE=0;
		for(my $k=1; $k<$a[8]; $k++) {
			my $Jid=$a[2]."_".$a[3]."_".$End[$k-1]."_".$Start[$k];
			my $Jstart=$End[$k-1];
			my $Jend=$Start[$k];
			if (exists $ModJunc{$Jid}) {
				my @b=split(/\_/,$ModJunc{$Jid});
				$Jstart=$b[-2];
				$Jend=$b[-1];
			}
			if (($force_illumina eq 0) and ((($NewJunc{$a[2]}{$a[3]}{$Jstart."\t".$Jend} >= $min_junc) or ((exists $KnownExonBorder{$a[2]}{$a[3]}{$Jstart}) and (exists $KnownExonBorder{$a[2]}{$a[3]}{$Jend}))) or (exists $SupportJunc{$a[2]}{$a[3]}{$Jstart."\t".$Jend}))) { }
			elsif (($force_illumina > 0) and (exists $KnownExonBorder{$a[2]}{$a[3]}{$Jstart}) and (exists $KnownExonBorder{$a[2]}{$a[3]}{$Jend}) and ($KnownExonBorder{$a[2]}{$a[3]}{$Jstart} >= $force_illumina) and ($KnownExonBorder{$a[2]}{$a[3]}{$Jend} >= $force_illumina)) {}
			else {
				$RARE++;
				if ($debug <= 20) {print OUTdebug join("\t","Rare",$id,$k,$Jid,$Jstart,$Jend),"\n";}
				last;
			}
		}
		if ($RARE > 0) { next; }
		my $start=100*(int($a[4]/100));
		my $end=100*(int($a[5]/100))+100;
		if (exists $Gene_loci{$a[2]}) {
			if (exists $Gene_loci{$a[2]}{$a[3]}) {
				my $flag=0;
				my %conflict;
				foreach my $pos(sort{$a <=> $b} keys %{$Gene_loci{$a[2]}{$a[3]}}) {
					if ($Gene_loci{$a[2]}{$a[3]}{$pos} < $start) {	}
					elsif ($end < $pos) { last; }
					else {
						$conflict{$pos}=$Gene_loci{$a[2]}{$a[3]}{$pos};
						$flag++;
					}
				}
				if ($flag eq 0) { $Gene_loci{$a[2]}{$a[3]}{$start}=$end; }
				else {
					my @Pos;
					push @Pos,$start;
					push @Pos,$end;
					foreach my $pos (keys %conflict) {
						push @Pos,$pos;
						push @Pos,$conflict{$pos};
						delete $Gene_loci{$a[2]}{$a[3]}{$pos};
					}
					my @POS=sort{$a <=> $b} @Pos;
					$Gene_loci{$a[2]}{$a[3]}{$POS[0]}=$POS[2*$flag+1];
				}
			}
			else { $Gene_loci{$a[2]}{$a[3]}{$start}=$end; }
		}
		else { $Gene_loci{$a[2]}{$a[3]}{$start}=$end; }
		if ($debug <= 50) {print OUTdebug join("\t",$a[2],$a[3],$start,$end,$Read{$a[0]}),"\n";}
	}
}

if($debug <= 20) {
	foreach my $id (sort keys %ModJunc) {
		my @a=split(/\_/,$id);
		my @b=split(/\_/,$ModJunc{$id});
		print OUTdebug "old\t",join("\t",@a),"\t",$Junc{$id},"\tnew\t",join("\t",@b),"\t",$Junc{$ModJunc{$id}},"\t",$NewJunc{$b[0]}{$b[1]}{$b[2]."\t".$b[3]},"\n";
	}
}
if ($debug <= 50){
	foreach my $chr (sort keys %Gene_loci) {
		foreach my $strand (sort keys %{$Gene_loci{$chr}}) {
			foreach my $start (sort{$a <=> $b} keys %{$Gene_loci{$chr}{$strand}}) {
				print OUTdebug "Gene_loci\t",join("_",$chr,$strand,$start,$Gene_loci{$chr}{$strand}{$start}),"\n";
			}
		}
	}
}

# assign multi-exon-reads and single-exon_reads_with_known_strandness to gene_loci
my %Loci;
foreach my $id (keys %Read) {
	if (exists $Single{$id}) { next; }
	else {
		my @a=split("\t",$Read{$id});
		foreach my $start (sort{$a <=> $b} keys %{$Gene_loci{$a[2]}{$a[3]}}) {
			my $end=$Gene_loci{$a[2]}{$a[3]}{$start};
			if (($start <= $a[4]) and ($a[5] <= $end)) {
				my $Lid=join("_",$a[2],$a[3],$start,$end);
				if (exists $Loci{$Lid}) { $Loci{$Lid}=$Loci{$Lid}."\t".$id;	}
				else { $Loci{$Lid}=$id; }
				last;
			} 
		}
	}
}
# filter gene_loci that consists only of rare reads
foreach my $Lid (sort keys %Loci) {
	my @L=split(/\_/,$Lid);
	my @Rid=split("\t",$Loci{$Lid});
	my $Nr=scalar(@Rid);
	my $RARE=0;
	for(my $i=0; $i<$Nr; $i++) {
		my @a=split("\t",$Read{$Rid[$i]});
		my @Start=split(/\,/,$a[9]);
		my @End=split(/\,/,$a[10]);
		for(my $k=1; $k<$a[8]; $k++) {
			my $id=$a[2]."_".$a[3]."_".$End[$k-1]."_".$Start[$k];
			my $start=$End[$k-1];
			my $end=$Start[$k];
			if (exists $ModJunc{$id}) {
				my @b=split(/\_/,$ModJunc{$id});
				$start=$b[-2];
				$end=$b[-1];
			}
			#if (($NewJunc{$L[0]}{$L[1]}{$start."\t".$end} >= $min_junc) or ((exists $KnownExonBorder{$L[0]}{$L[1]}{$start}) and (exists $KnownExonBorder{$L[0]}{$L[1]}{$end}))) { }
			if (($force_illumina eq 0) and ((($NewJunc{$L[0]}{$L[1]}{$start."\t".$end} >= $min_junc) or ((exists $KnownExonBorder{$L[0]}{$L[1]}{$start}) and (exists $KnownExonBorder{$L[0]}{$L[1]}{$end}))) or (exists $SupportJunc{$L[0]}{$L[1]}{$start."\t".$end}))) { }
			elsif (($force_illumina > 0) and (exists $KnownExonBorder{$L[0]}{$L[1]}{$start}) and (exists $KnownExonBorder{$L[0]}{$L[1]}{$end}) and ($KnownExonBorder{$L[0]}{$L[1]}{$start} >= $force_illumina) and ($KnownExonBorder{$L[0]}{$L[1]}{$end} >= $force_illumina)) {}
			else {
				$RareRead{$Rid[$i]}=1;
				$RARE++;
				last;
			}
		}
	}
	if ($RARE eq $Nr) {
		# this loci consists only of rare reads
		delete $Loci{$Lid};
	}
}

# assign single-exon-read_without_strandness to known_gene_loci
foreach my $id (keys %Read) {
	if (exists $Single{$id}) {
		my @a=split("\t",$Read{$id});
		my $flagp="na";
		my $flagm="na";
		my $strand="+";
		# check if overlap with known_gene_loci
		foreach my $start (sort{$a <=> $b} keys %{$Gene_loci{$a[2]}{$strand}}) {
			my $end=$Gene_loci{$a[2]}{$strand}{$start};
			if (($a[5] < ($start - $fuzzy)) or ($a[4] > ($end + $fuzzy))) { }
			else {
			#if (($start <= $a[4]) and ($a[5] <= $end)) {
				my $Lid=join("_",$a[2],$strand,$start,$end);
				$flagp=$Lid;
				last;
			} 
		}
		if ($debug <= 50){ print OUTdebug join("\t","LineF504",$id,$a[2],$a[3],$strand,$flagp),"\n"; }
		$strand="-";
		foreach my $start (sort{$a <=> $b} keys %{$Gene_loci{$a[2]}{$strand}}) {
			my $end=$Gene_loci{$a[2]}{$strand}{$start};
			if (($a[5] < ($start - $fuzzy)) or ($a[4] > ($end + $fuzzy))) { }
			else {
			#if (($start <= $a[4]) and ($a[5] <= $end)) {
				my $Lid=join("_",$a[2],$strand,$start,$end);
				$flagm=$Lid;
				last;
			} 
		}
		if ($debug <= 50){ print OUTdebug join("\t","LineF516",$id,$a[2],$a[3],$strand,$flagm),"\n"; }
		if (($flagp ne "na") and ($flagm ne "na")) {
			$strand="na";
			my $flag=0;
			foreach my $start (sort{$a <=> $b} keys %{$Gene_loci{$a[2]}{$strand}}) {
				my $end=$Gene_loci{$a[2]}{$strand}{$start};
				if (($start <= $a[4]) and ($a[5] <= $end)) {
					$flag++;
					my $Lid=join("_",$a[2],$strand,$start,$end);
					if (exists $Loci{$Lid}) { $Loci{$Lid}=$Loci{$Lid}."\t".$id;	}
					else { $Loci{$Lid}=$id; }
					last;
				} 
			}
			if ($flag eq 0) {
				# this read doesn't fall within any known_gene_loci, check if this read overlap with any neutral_gene_loci
				foreach my $start (sort{$a <=> $b} keys %{$Gene_loci{$a[2]}{$strand}}) {
					my $end=$Gene_loci{$a[2]}{$strand}{$start};
					if (($a[5] < ($start - $fuzzy)) or ($a[4] > ($end + $fuzzy))) { }
					else {
					#if (($start <= $a[4]) and ($a[5] <= $end)) {
						$flag++;
						my $Lid=join("_",$a[2],$strand,$start,$end);
						if (exists $Loci{$Lid}) { $Loci{$Lid}=$Loci{$Lid}."\t".$id;	}
						else { $Loci{$Lid}=$id; }
						last;
					} 
				}
				if ($flag eq 0) {
					# this read doesn't fall within any known_gene_loci, create another neutral_gene_locus
					my $start=100*(int($a[4]/100));
					my $end=100*(int($a[5]/100))+100;
					my $Lid=join("_",$a[2],$strand,$start,$end);
					$Loci{$Lid}=$id;
					$Gene_loci{$a[2]}{$strand}{$start}=$end;
				}
			}
		}
		elsif (($flagp eq "na") and ($flagm eq "na")) {
			# check if overlap with neutral_gene_loci
			$strand="na";
			my $flag=0;
			foreach my $start (sort{$a <=> $b} keys %{$Gene_loci{$a[2]}{$strand}}) {
				my $end=$Gene_loci{$a[2]}{$strand}{$start};
				if (($a[5] < ($start - $fuzzy)) or ($a[4] > ($end + $fuzzy))) { }
				else {
				#if (($start <= $a[4]) and ($a[5] <= $end)) {
					$flag++;
					my $Lid=join("_",$a[2],$strand,$start,$end);
					if (exists $Loci{$Lid}) { $Loci{$Lid}=$Loci{$Lid}."\t".$id;	}
					else { $Loci{$Lid}=$id; }
					last;
				} 
			}
			if ($flag eq 0) {
				# this read doesn't fall within any known_gene_loci, create another neutral_gene_locus
				my $start=100*(int($a[4]/100));
				my $end=100*(int($a[5]/100))+100;
				my $Lid=join("_",$a[2],$strand,$start,$end);
				$Loci{$Lid}=$id;
				$Gene_loci{$a[2]}{$strand}{$start}=$end;
			}
		}
		elsif ($flagp ne "na") {
			if (exists $Loci{$flagp}) { $Loci{$flagp}=$Loci{$flagp}."\t".$id;	}
			else { $Loci{$flagp}=$id; }
		}
		elsif ($flagm ne "na") {
			if (exists $Loci{$flagm}) { $Loci{$flagm}=$Loci{$flagm}."\t".$id;	}
			else { $Loci{$flagm}=$id; }
		}
		if ($debug <= 50){ print OUTdebug join("\t","LineF568",$Read{$id},$flagp,$flagm,$strand),"\n";	}
		
	}
}



if ($debug <= 30) {
	foreach my $Lid (sort keys %Loci) {
		my @a=split("\t",$Loci{$Lid});
		my $info=$Lid."\t".scalar(@a);
		for(@a) {$info=$info."\t".$_};
		print OUTdebug $info,"\n";
	}
}


# foreach gene_loci, cluster reads according to junc
my $ClusterCnt=0;

foreach my $Lid (sort keys %Loci) {
	my $yesprint=0;
	$ClusterCnt++;
	if ($debug <= 50) { print OUTdebug $ClusterCnt,"\t",$Lid,"\t",$Loci{$Lid},"\n";	}
	my @L=split(/\_/,$Lid);
	my @Rid=split("\t",$Loci{$Lid});
	my $Nr=scalar(@Rid);
	if ($Nr eq 1) {
		# output this single track
		my @a=split("\t",$Read{$Rid[0]});
		$a[0]="GL.".$ClusterCnt;
		print OUT6 join("\t",@a),"\n";
		next;
	}
	if ($L[1] eq "+") {
		# +strand, 5'UTR on the left-most
		# make a survey of all junctions within this loci
		my %Sjunc;
		my %Five;
		my %Three;
		my %Readjunc;
		for(my $i=0; $i<$Nr; $i++) {
			my @a=split("\t",$Read{$Rid[$i]});
			if (exists $RareRead{$Rid[$i]}) {next;};
			$Five{$Rid[$i]}=$a[4];
			$Three{$Rid[$i]}=$a[5];
			my @Start=split(/\,/,$a[9]);
			my @End=split(/\,/,$a[10]);
			for(my $j=1; $j<$a[8]; $j++) {
				my $jid=$L[0]."_".$L[1]."_".$End[$j-1]."_".$Start[$j];
				if (exists $ModJunc{$jid}) { my $tmp=$ModJunc{$jid}; $jid=$tmp; }
				my @t=split(/\_/,$jid);
				$Sjunc{$t[2]}{$t[3]}=1;
				$Readjunc{$Rid[$i]}{$jid}=1;
			}
		}
		# make left-to-right full-pattern
		my @SJunc;
		my $cnt=0;
		foreach my $start (sort{$a <=> $b} keys %Sjunc) {
			foreach my $end (sort{$a <=> $b} keys %{$Sjunc{$start}} ) {
				# filter rare junctions which are likely to be wrong
				if (($force_illumina eq 0) and ((($NewJunc{$L[0]}{$L[1]}{$start."\t".$end} >= $min_junc) or ((exists $KnownExonBorder{$L[0]}{$L[1]}{$start}) and (exists $KnownExonBorder{$L[0]}{$L[1]}{$end}))) or (exists $SupportJunc{$L[0]}{$L[1]}{$start."\t".$end}))) {
					my $jid=join("_",$L[0],$L[1],$start,$end);
					$SJunc[$cnt]=$jid;
					$cnt++;
				}
				elsif (($force_illumina > 0) and (exists $KnownExonBorder{$L[0]}{$L[1]}{$start}) and (exists $KnownExonBorder{$L[0]}{$L[1]}{$end}) and ($KnownExonBorder{$L[0]}{$L[1]}{$start} >= $force_illumina) and ($KnownExonBorder{$L[0]}{$L[1]}{$end} >= $force_illumina)) {
					my $jid=join("_",$L[0],$L[1],$start,$end);
					$SJunc[$cnt]=$jid;
					$cnt++;
				}
				#if (($NewJunc{$L[0]}{$L[1]}{$start."\t".$end} >= $min_junc) or ((exists $KnownExonBorder{$L[0]}{$L[1]}{$start}) and (exists $KnownExonBorder{$L[0]}{$L[1]}{$end}))) {
				#	my $jid=join("_",$L[0],$L[1],$start,$end);
				#	$SJunc[$cnt]=$jid;
				#	$cnt++;
				#}
				else {
					foreach my $rid (keys %Readjunc) {
						my $jid=join("_",$L[0],$L[1],$start,$end);
						if (exists $Readjunc{$rid}{$jid}) {
							$RareRead{$rid}=1;
						}
					}
				}
			}
		}
		# scan for each Read and make pattern
		my %pattern;
		my %store;
		my $mscore=-1;
		my $blank=0;
		for(my $i=1; $i<$cnt; $i++) {$blank=$blank."0";}
		for(my $i=0; $i<$Nr; $i++) {
			if (exists $RareRead{$Rid[$i]}) { next; }
			my $p="";
			my $start=0;		# Read's start wrt pattern
			my $end=$cnt-1;		# Read's end   wrt pattern
			my $score=0;
			for(my $j=0; $j<$cnt; $j++) {
				if (exists $Readjunc{$Rid[$i]}{$SJunc[$j]}) {
					my @ta=split(/\_/,$SJunc[$j]);
					if (($Five{$Rid[$i]} + $fuzzy)  > $ta[3]) { $start=$j+1; }
					if (($Three{$Rid[$i]} < $ta[2]) and ($end eq ($cnt-1))) { $end=$j-1; }
					$p=$p."1";
					$score+=2**$j;
				}
				else {
					$p=$p."0";
					my @ta=split(/\_/,$SJunc[$j]);
					if (($Five{$Rid[$i]} + $fuzzy)  > $ta[3]) { $start=$j+1; }
					if (($Three{$Rid[$i]} < $ta[2]) and ($end eq ($cnt-1))) { $end=$j-1; }
				}
			}
			if ($start >= $end) {
				# still report his
				## ignore single exon reads if the geneLoci has >=2 junctions
				#if ($cnt > 1) {
				#	next;
				#}
				# $end=$start;
			}
			$pattern{$end-$start+1}{$score}{$Len{$Rid[$i]}}{$Rid[$i]}=join("\t",$p,$start,$end,$score);
			$store{$Rid[$i]}=join("\t",$p,$start,$end,$score);
			if ($score > $mscore) {$mscore = $score;}
		}
		# scan for each pattern and cluster
		my $NrCluster=0;
		my @MyCluster;
		my %myCluster;
		foreach my $span (sort{$b <=> $a} keys %pattern) {
			foreach my $score (sort{$b <=> $a} keys %{$pattern{$span}}) {
				foreach my $len (sort{$b <=> $a} keys %{$pattern{$span}{$score}}) {
					foreach my $rid (sort keys %{$pattern{$span}{$score}{$len}}) {
						if ($debug <= 60) {print "this == ",$rid,"\tNrCluster = ",$NrCluster,"\t",$span,"\t",$score,"\t",$len,"\t",$pattern{$span}{$score}{$len}{$rid},"\n";}
						my @a=split("\t",$pattern{$span}{$score}{$len}{$rid});
						# check for truncated reads
						#if (($mscore > 1) and ($a[0] eq $blank)) {
						if (($a[0] eq "") or ($a[0] eq $blank)) {
							$Truncate{$rid}=1;
							next;
						}
						if ($NrCluster eq 0) {
							$NrCluster=1;
							$MyCluster[$NrCluster]=join("\t",$rid,$a[0],$a[1],$a[2],$a[3]);
							$myCluster{$rid}=$rid;
						}
						else {
							my $flag=-1;
							for (my $k=1; $k<=$NrCluster; $k++) {
								my @b=split("\t",$MyCluster[$k]);
								my $stringa=substr($a[0],$a[1],($a[2]-$a[1]+1));
								my $stringb=substr($b[1],$a[1],($a[2]-$a[1]+1));
								if ($stringa eq $stringb) {
									# @a is a subset of $MyCluster[$k]
									if ($flag eq -1) { $flag=$k; }
									else {$flag=$flag."\t".$k;}
								}
								if ($debug <= 60) {
									print $rid,"\t",join("\t",@a),"\tstringa == ",$stringa,"\t",$Five{$rid},"\t",$Three{$rid},"\n";
									print $MyCluster[$k],"\tstringb == ",$stringb,"\t",$Five{$b[0]},"\t",$Three{$b[0]},"\n";
								}
							}
							if ($debug <= 60) {print "flag == ",$flag,"\n";}
							if ($flag eq -1) {
								if (!exists $Truncate{$rid}) {
									$NrCluster++;
									$MyCluster[$NrCluster]=join("\t",$rid,$a[0],$a[1],$a[2],$a[3]);
									$myCluster{$rid}=$rid;
								}
							}
							else {
								my $same=-1;
								my @flaga=split("\t",$flag);
								my $Nrfa=scalar(@flaga);
								for(my $k=0; ($k<$Nrfa)and($same eq -1); $k++) {
									my @b=split("\t",$MyCluster[$flaga[$k]]);
									# if 3'end is near
									if (abs($Three{$rid} - $Three{$b[0]}) < $min_3UTR_diff) {
										# check if this 5'end overlap with Cage
										my @ra=split("\t",$Read{$rid});
										my $rbin=int($ra[4]/1000);
										my $altStart=0;
										if (($knowncage ne "no") and (exists $KnownCage{$L[0]}{$L[1]}{$rbin})) {
											foreach my $CageStart (keys %{$KnownCage{$L[0]}{$L[1]}{$rbin}}) {
												if (abs($CageStart - $Five{$rid}) < 100) {
													$altStart=$CageStart;
													last;
												}
											}
										}
										if (($altStart eq 0) and ($knowncage ne "no") and (exists $KnownCage{$L[0]}{$L[1]}{$rbin+1})) {
											foreach my $CageStart (keys %{$KnownCage{$L[0]}{$L[1]}{$rbin+1}}) {
												if (abs($CageStart - $Five{$rid}) < 100) {
													$altStart=$CageStart;
													last;
												}
											}
										}
										if (($altStart eq 0) and ($knowncage ne "no") and (exists $KnownCage{$L[0]}{$L[1]}{$rbin-1})) {
											foreach my $CageStart (keys %{$KnownCage{$L[0]}{$L[1]}{$rbin-1}}) {
												if (abs($CageStart - $Five{$rid}) < 100) {
													$altStart=$CageStart;
													last;
												}
											}
										}
										if ($altStart ne 0) {
											# Overlap with known Start
										}
										else {
											# check if 5'end is within exon
											#my @c=split("\t",$Read{$b[0]});
											#my @CStart=split(/\,/,$c[9]);
											#my @CEnd=split(/\,/,$c[10]);
											#my $within=-1;
											#for(my $kk=0; ($kk<$c[8])and($within eq -1); $kk++) {
											#	if (($Five{$rid} > ($CStart[$kk] - $fuzzy)) and ($Five{$rid} < ($CEnd[$kk]+$fuzzy))) {
											#		$within=$kk;
											#	}
											#}
											#if ($within ne -1) {
											#	# 5'end within known exons
											#	$same=$flaga[$k];
											#}
											# check if 5' end is within any exon
											my $within=-1;
											foreach my $kid (keys %store) {
												my @c=split("\t",$Read{$kid});
												my @CStart=split(/\,/,$c[9]);
												my @CEnd=split(/\,/,$c[10]);
												for(my $kk=0; ($kk<$c[8])and($within eq -1); $kk++) {
													if (($Five{$rid} > ($CStart[$kk] - $fuzzy)) and ($Five{$rid} < ($CEnd[$kk]+$fuzzy))) {
														$within=$kk;
													}
												}
											}
											if ($within ne -1) {
												# 5'end within known exons
												$same=$flaga[$k];
											}
										}
									}
									# if 3'end is far away, possible new isofrom
									else {
										# 
									}
								}
								if ($same eq -1) {
									$NrCluster++;
									$MyCluster[$NrCluster]=join("\t",$rid,$a[0],$a[1],$a[2],$a[3]);
									$myCluster{$rid}=$rid;
								}
								else {
									my @tmp=split("\t",$MyCluster[$same]);
									$myCluster{$rid}=$tmp[0];
								}
							}
						}
					}
				}
			}
		}
		# report
		for (my $i=1; $i<=$NrCluster; $i++) {
			my @a=split("\t",$MyCluster[$i]);
			my @b=split("\t",$Read{$a[0]});
			$b[0]="GL.".$ClusterCnt;
			# fix junctions
			my @Start=split(/\,/,$b[9]);
			my @End=split(/\,/,$b[10]);
			for(my $j=1; $j<$b[8]; $j++) {
				my $jid=$L[0]."_".$L[1]."_".$End[$j-1]."_".$Start[$j];
				if (exists $ModJunc{$jid}) {
					my @c=split(/\_/,$ModJunc{$jid});
					$End[$j-1]=$c[2];
					$Start[$j]=$c[3];
				}
			}
			# fill in $b[11] for all the sub-isofroms
			my $tail=$b[1];
			if ($debug <= 40) {
				print OUTdebug join("\t",$b[0],$MyCluster[$i],$Five{$b[1]},$Three{$b[1]}),"\n";
			}			
			foreach my $trid (keys %myCluster) {
				if (($trid ne $b[1]) and ($myCluster{$trid} eq $b[1])) {
					$tail=$tail.",".$trid;
					if ($Start[0] > $Five{$trid}) { $Start[0]=$Five{$trid}; }
					if ($End[-1] < $Three{$trid}) { $End[-1]=$Three{$trid}; }
					if ($debug <= 40) {
						print OUTdebug join("\t",$b[0].".s",$trid,$store{$trid},$Five{$trid},$Three{$trid}),"\n";
					}
				}
			}
			$b[4]=$Start[0];
			$b[5]=$End[-1];
			$b[6]=$b[4];
			$b[7]=$b[4];
			# sanity check for exon structures
			my $sany=0;
			if ($End[0] <= $Start[0]) { $sany++; }
			for(my $i=1; $i<$b[8]; $i++) {
				if ($Start[$i] < $End[$i-1]) { $sany++; }
				if ($End[$i] < $Start[$i]) { $sany++; }
			}
			$b[9]=join(",",@Start);
			$b[10]=join(",",@End);
			push @b, $tail;
			if ($sany > 0) {
				print OUT10 join("\t",@b),"\n";
			}
			else {
				print OUT1 join("\t",@b),"\n";
				$yesprint++;
			}
		}
	}
	elsif($L[1] eq "-") {
		# -strand, 5'UTR on the right-most
		# make a survey of all junctions within this loci
		my %Sjunc;
		my %Five;
		my %Three;
		my %Readjunc;
		for(my $i=0; $i<$Nr; $i++) {
			my @a=split("\t",$Read{$Rid[$i]});
			if (exists $RareRead{$Rid[$i]}) {next;};
			$Five{$Rid[$i]}=$a[5];
			$Three{$Rid[$i]}=$a[4];
			my @Start=split(/\,/,$a[9]);
			my @End=split(/\,/,$a[10]);
			for(my $j=1; $j<$a[8]; $j++) {
				my $jid=$L[0]."_".$L[1]."_".$End[$j-1]."_".$Start[$j];
				if (exists $ModJunc{$jid}) { my $tmp=$ModJunc{$jid}; $jid=$tmp; }
				my @t=split(/\_/,$jid);
				$Sjunc{$t[2]}{$t[3]}=1;
				$Readjunc{$Rid[$i]}{$jid}=1;
			}
		}
		# make left-to-right full-pattern
		my @SJunc;
		my $cnt=0;
		foreach my $start (sort{$a <=> $b} keys %Sjunc) {
			foreach my $end (sort{$a <=> $b} keys %{$Sjunc{$start}} ) {
				# filter rare junctions which are likely to be wrong
				if (($force_illumina eq 0) and ((($NewJunc{$L[0]}{$L[1]}{$start."\t".$end} >= $min_junc) or ((exists $KnownExonBorder{$L[0]}{$L[1]}{$start}) and (exists $KnownExonBorder{$L[0]}{$L[1]}{$end}))) or (exists $SupportJunc{$L[0]}{$L[1]}{$start."\t".$end}))) {
					my $jid=join("_",$L[0],$L[1],$start,$end);
					$SJunc[$cnt]=$jid;
					$cnt++;
				}
				elsif (($force_illumina > 0) and (exists $KnownExonBorder{$L[0]}{$L[1]}{$start}) and (exists $KnownExonBorder{$L[0]}{$L[1]}{$end}) and ($KnownExonBorder{$L[0]}{$L[1]}{$start} >= $force_illumina) and ($KnownExonBorder{$L[0]}{$L[1]}{$end} >= $force_illumina)) {
					my $jid=join("_",$L[0],$L[1],$start,$end);
					$SJunc[$cnt]=$jid;
					$cnt++;
				}
				#if (($NewJunc{$L[0]}{$L[1]}{$start."\t".$end} >= $min_junc) or ((exists $KnownExonBorder{$L[0]}{$L[1]}{$start}) and (exists $KnownExonBorder{$L[0]}{$L[1]}{$end}))) {
				#	my $jid=join("_",$L[0],$L[1],$start,$end);
				#	$SJunc[$cnt]=$jid;
				#	$cnt++;
				#}
				else {
					foreach my $rid (keys %Readjunc) {
						my $jid=join("_",$L[0],$L[1],$start,$end);
						if (exists $Readjunc{$rid}{$jid}) {
							$RareRead{$rid}=1;
						}
					}
				}
			}
		}
		# scan for each Read and make pattern
		my %pattern;
		my %store;
		my $mscore=-1;
		my $blank=0;
		for(my $i=1; $i<$cnt; $i++) {$blank=$blank."0";}
		for(my $i=0; $i<$Nr; $i++) {
			if (exists $RareRead{$Rid[$i]}) { next; }
			my $p="";
			my $start=$cnt-1;		# Read's start wrt pattern
			my $end=0;				# Read's end   wrt pattern
			my $score=0;
			for(my $j=0; $j<$cnt; $j++) {
				if (exists $Readjunc{$Rid[$i]}{$SJunc[$j]}) {
					my @ta=split(/\_/,$SJunc[$j]);
					if ($Three{$Rid[$i]} > $ta[3]) { $end=$j+1; }
					if ((($Five{$Rid[$i]} - $fuzzy) < $ta[2]) and ($start eq ($cnt-1))) { $start=$j-1; }
					$p=$p."1";
					$score+=2**($cnt-$j);
				}
				else {
					$p=$p."0";
					my @ta=split(/\_/,$SJunc[$j]);
					if ($Three{$Rid[$i]} > $ta[3]) { $end=$j+1; }
					if ((($Five{$Rid[$i]} - $fuzzy) < $ta[2]) and ($start eq ($cnt-1))) { $start=$j-1; }
				}
			}
			if ($end >= $start) {
				# still report his
				## ignore single exon reads if the geneLoci has >=2 junctions
				#if ($cnt > 1) {
				#	next;
				#}
				#$end=$start;
			}
			$pattern{$start-$end+1}{$score}{$Len{$Rid[$i]}}{$Rid[$i]}=join("\t",$p,$start,$end,$score);
			$store{$Rid[$i]}=join("\t",$p,$start,$end,$score);
			if ($score > $mscore) {$mscore = $score;}
		}
		# scan for each pattern and cluster
		my $NrCluster=0;
		my @MyCluster;
		my %myCluster;
		foreach my $span (sort{$b <=> $a} keys %pattern) {
			foreach my $score (sort{$b <=> $a} keys %{$pattern{$span}}) {
				foreach my $len (sort{$b <=> $a} keys %{$pattern{$span}{$score}}) {
					foreach my $rid (sort keys %{$pattern{$span}{$score}{$len}}) {
						if ($debug <= 60) {print "this == ",$rid,"\tNrCluster = ",$NrCluster,"\t",$span,"\t",$score,"\t",$len,"\t",$pattern{$span}{$score}{$len}{$rid},"\n";}
						my @a=split("\t",$pattern{$span}{$score}{$len}{$rid});
						# check for truncated reads
						# if (($mscore > 1) and ($a[0] eq $blank)) {
						if (($a[0] eq "") or ($a[0] eq $blank)) {
							$Truncate{$rid}=1;
							next;
						}
						if ($NrCluster eq 0) {
							$NrCluster=1;
							$MyCluster[$NrCluster]=join("\t",$rid,$a[0],$a[1],$a[2],$a[3]);
							$myCluster{$rid}=$rid;
						}
						else {
							my $flag=-1;
							for (my $k=1; $k<=$NrCluster; $k++) {
								my @b=split("\t",$MyCluster[$k]);
								my $stringa=substr($a[0],$a[2],($a[1]-$a[2]+1));
								my $stringb=substr($b[1],$a[2],($a[1]-$a[2]+1));
								if ($stringa eq $stringb) {
									# if they are both single-exon-reads, most likely to be truncated 3'UTR
									if (($stringa eq "0") or ($stringa eq "")) {
										$Truncate{$rid}=1;
										# print OUT4 join("\t",$Lid,join("\t",@a),$Read{$rid}),"\n";
									}
									# @a is a subset of $MyCluster[$k]
									else {
										if ($flag eq -1) { $flag=$k; }
										else {$flag=$flag."\t".$k;}
									}
								}
								if ($debug <= 60) {
									print $rid,"\t",join("\t",@a),"\tstringa == ",$stringa,"\t",$Five{$rid},"\t",$Three{$rid},"\n";
									print $MyCluster[$k],"\tstringb == ",$stringb,"\t",$Five{$b[0]},"\t",$Three{$b[0]},"\n";
								}
							}
							if ($debug <= 60) {print "flag == ",$flag,"\n";}
							if ($flag eq -1) {
								if (!exists $Truncate{$rid}) {
									$NrCluster++;
									$MyCluster[$NrCluster]=join("\t",$rid,$a[0],$a[1],$a[2],$a[3]);
									$myCluster{$rid}=$rid;
								}
							}
							else {
								my $same=-1;
								my @flaga=split("\t",$flag);
								my $Nrfa=scalar(@flaga);
								for(my $k=0; ($k<$Nrfa)and($same eq -1); $k++) {
									my @b=split("\t",$MyCluster[$flaga[$k]]);
									# if 3'end is near
									if ($debug <= 61) {
										print join("\t",$same,$rid,join("\t",@a),$Five{$rid},$Three{$rid}),"\n";
										print join("\t",$same,$MyCluster[$flaga[$k]],$Five{$b[0]},$Three{$b[0]}),"\n";
									}
									if (abs($Three{$rid} - $Three{$b[0]}) < $min_3UTR_diff) {
										# check if this 5'end overlap with Cage
										my @ra=split("\t",$Read{$rid});
										my $rbin=int($ra[4]/1000);
										my $altStart=0;
										if (($knowncage ne "no") and (exists $KnownCage{$L[0]}{$L[1]}{$rbin})) {
											foreach my $CageStart (keys %{$KnownCage{$L[0]}{$L[1]}{$rbin}}) {
												if (abs($CageStart - $Five{$rid}) < 100) {
													$altStart=$CageStart;
													last;
												}
											}
										}
										if (($altStart eq 0) and ($knowncage ne "no") and (exists $KnownCage{$L[0]}{$L[1]}{$rbin+1})) {
											foreach my $CageStart (keys %{$KnownCage{$L[0]}{$L[1]}{$rbin+1}}) {
												if (abs($CageStart - $Five{$rid}) < 100) {
													$altStart=$CageStart;
													last;
												}
											}
										}
										if (($altStart eq 0) and ($knowncage ne "no") and (exists $KnownCage{$L[0]}{$L[1]}{$rbin-1})) {
											foreach my $CageStart (keys %{$KnownCage{$L[0]}{$L[1]}{$rbin-1}}) {
												if (abs($CageStart - $Five{$rid}) < 100) {
													$altStart=$CageStart;
													last;
												}
											}
										}
										if ($altStart ne 0) {
											# Overlap with known Start
										}
										else {
											# check if 5' end is within any exon
											my $within=-1;
											foreach my $kid (keys %store) {
												my @c=split("\t",$Read{$kid});
												my @CStart=split(/\,/,$c[9]);
												my @CEnd=split(/\,/,$c[10]);
												for(my $kk=0; ($kk<$c[8])and($within eq -1); $kk++) {
													if (($Five{$rid} > ($CStart[$kk] - $fuzzy)) and ($Five{$rid} < ($CEnd[$kk]+$fuzzy))) {
														$within=$kk;
													}
												}
											}
											if ($within ne -1) {
												# 5'end within known exons
												$same=$flaga[$k];
											}
										}
									}
									# if 3'end is far away, possible new isofrom
									else {
										# 
									}
								}
								if ($same eq -1) {
									$NrCluster++;
									$MyCluster[$NrCluster]=join("\t",$rid,$a[0],$a[1],$a[2],$a[3]);
									$myCluster{$rid}=$rid;
								}
								else {
									my @tmp=split("\t",$MyCluster[$same]);
									$myCluster{$rid}=$tmp[0];
								}
							}
						}
					}
				}
			}
		}
		# report
		for (my $i=1; $i<=$NrCluster; $i++) {
			my @a=split("\t",$MyCluster[$i]);
			my @b=split("\t",$Read{$a[0]});
			$b[0]="GL.".$ClusterCnt;
			# fix junctions
			my @Start=split(/\,/,$b[9]);
			my @End=split(/\,/,$b[10]);
			for(my $j=1; $j<$b[8]; $j++) {
				my $jid=$L[0]."_".$L[1]."_".$End[$j-1]."_".$Start[$j];
				if (exists $ModJunc{$jid}) {
					my @c=split(/\_/,$ModJunc{$jid});
					$End[$j-1]=$c[2];
					$Start[$j]=$c[3];
				}
			}			
			# fill in $b[11] for all the sub-isofroms
			my $tail=$b[1];
			if ($debug <= 40) {
				print OUTdebug join("\t",$b[0],$MyCluster[$i],$Five{$b[1]},$Three{$b[1]}),"\n";
			}			
			foreach my $trid (keys %myCluster) {
				if (($trid ne $b[1]) and ($myCluster{$trid} eq $b[1])) {
					$tail=$tail.",".$trid;
					if ($Start[0] > $Five{$trid}) { $Start[0]=$Five{$trid}; }
					if ($End[-1] < $Three{$trid}) { $End[-1]=$Three{$trid}; }
					if ($debug <= 40) {
						print OUTdebug join("\t",$b[0].".s",$trid,$store{$trid},$Five{$trid},$Three{$trid}),"\n";
					}
				}
			}
			$b[4]=$Start[0];
			$b[5]=$End[-1];
			$b[6]=$b[4];
			$b[7]=$b[4];
			# sanity check for exon structures
			my $sany=0;
			if ($End[0] <= $Start[0]) { $sany++; }
			for(my $i=1; $i<$b[8]; $i++) {
				if ($Start[$i] < $End[$i-1]) { $sany++; }
				if ($End[$i] < $Start[$i]) { $sany++; }
			}
			$b[9]=join(",",@Start);
			$b[10]=join(",",@End);
			push @b, $tail;
			if ($sany > 0) {
				print OUT10 join("\t",@b),"\n";
			}
			else {
				print OUT1 join("\t",@b),"\n";
				$yesprint++;
			}
		}
	}
	else {
		# single-exon-read without any strandness information
		my %Five;
		my %FiveRid;
		my %Three;
		my %ThreeRid;
		my $left_most=-1;
		my $right_most=-1;
		for(my $i=0; $i<$Nr; $i++) {
			my @a=split("\t",$Read{$Rid[$i]});
			if ($left_most eq -1) { $left_most=$a[4]; }
			else { if($left_most > $a[4]){ $left_most=$a[4]; } }
			if ($right_most eq -1) { $right_most=$a[5]; }
			else { if($right_most < $a[5]){ $right_most=$a[5]; } }
			my $bin=int($a[4]/50);
			$Five{$bin}++;
			if (exists $FiveRid{$bin}) { $FiveRid{$bin}=$FiveRid{$bin}."\t".$Rid[$i]; }
			else { $FiveRid{$bin}=$Rid[$i]; }
			$bin=int($a[5]/50);
			$Three{$bin}++;
			if (exists $ThreeRid{$bin}) { $ThreeRid{$bin}=$ThreeRid{$bin}."\t".$Rid[$i]; }
			else { $ThreeRid{$bin}=$Rid[$i]; }
		}
		my $Fivem=0;
		foreach my $pos (sort{$Five{$b} <=> $Five{$a}} keys %Five) {
			$Fivem=$pos;
			last;
		}
		my $Threem=0;
		foreach my $pos (sort{$Three{$b} <=> $Three{$a}} keys %Three) {
			$Threem=$pos;
			last;
		}
		my @a=split("\t",$FiveRid{$Fivem});
		my $most="";
		for(@a) {
			my $tid=$_;
			if ($ThreeRid{$Threem}=~m/$tid/) { $most=$tid; last; }
		}
		if ($most eq "") {
			my @a=split("\t",$Read{$Rid[0]});
			$a[0]="GL.".$ClusterCnt;
			$a[4]=$left_most;
			$a[5]=$right_most;
			if($a[3] eq "+") {$a[6] = $left_most; $a[7] = $left_most;}
			if($a[3] eq "-") {$a[6] = $right_most; $a[7] = $right_most;}
			$a[9]=$left_most;
			$a[10]=$right_most;
			my $tail=join("\,",@Rid);
			#my $tail=$Rid[0];
			#for(my $i=1; $i<$Nr; $i++) { $tail=$tail.",".$Rid[$i]; }
			push @a, $tail;
			print OUT1 join("\t",@a),"\n";
		}
		else {
			my @a=split("\t",$Read{$most});
			$a[0]="GL.".$ClusterCnt;
			if($a[3] eq "+") {$a[6] = $a[4]; $a[7] = $a[4];}
			if($a[3] eq "-") {$a[6] = $a[5]; $a[7] = $a[5];}
			$a[4]=$left_most;
			$a[5]=$right_most;
			$a[9]=$left_most;
			$a[10]=$right_most;
			my $tail=join("\,",@Rid);
			#my $tail=$Rid[0];
			#for(my $i=1; $i<$Nr; $i++) { $tail=$tail.",".$Rid[$i]; }
			push @a, $tail;
			print OUT1 join("\t",@a),"\n";
		}
		$yesprint++;
	}
	if ($yesprint > 0) {
		my @a=split(/\_/,$Lid);
		# print OUT0 join("\t",$a[0],$a[2],$a[3],$Lid,$yesprint,$a[1]),"\n";
		print OUT0 join("\t",$a[0],$a[2],$a[3],"GL.".$ClusterCnt,$yesprint,$a[1],$Nr),"\n";
	}
}
foreach my $id (sort keys %Truncate) { print OUT4 $Read{$id},"\n"; }
foreach my $id (sort keys %RareRead) { print OUT5 $Read{$id},"\n"; }