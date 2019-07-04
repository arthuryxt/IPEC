#!/usr/bin/perl -w
die "Usage: $0  \"input_fa\"   \"3p_adapter\"   \"5p_adapter\"   \"\(optional\)error_rate : between 0 and 1; 0.1 by default\"   \"\(optional\) minLen : 200 by default \"   \"\(optional\) verbose = 1\"   \"\(optional\) remakeID = 1\" " if (@ARGV < 3);
my $filein1=$ARGV[0];
my $adapter_3p=$ARGV[1];    # AAGCAGTGGTATCAACGCAGAGTACTTTTT
my $adapter_5p=$ARGV[2];    # AAGCAGTGGTATCAACGCAGAGTACATGGG
my $ErrorRate=0.1;	    # 0.3 for pacbio
if (scalar(@ARGV) > 3) {$ErrorRate=$ARGV[3];}
my $minLen=300;
if (scalar(@ARGV) > 4) {$minLen=$ARGV[4];}
my $verbose=0;
if (scalar(@ARGV) > 5) {$verbose=$ARGV[5];}
my $remakeID=0;
if (scalar(@ARGV) > 6) {$remakeID=$ARGV[6];}


### by arthur 2013-09-20
### 5p_adapter --- [ 5p_cNDA_primer --- plus_strand_read_1 --- 3p_cDNA_primer ] --- 3p_adapter ###
### AAGCAGTGGTATCAACGCAGAGTACTTTTT --- TTTNNNNNN --- CCCATGTACTCTGCGTTGATACCACTGCTT
### AAGCAGTGGTATCAACGCAGAGTACATGGG --- NNNNNNAAA --- AAAAAGTACTCTGCGTTGATACCACTGCTT

# from Clontech
# my $adapter_5p="AAGCAGTGGTATCAACGCAGAGTACATGGG";	# F1 or F0
# my $adapter_3p="GTACTCTGCGTTGATACCACTGCTT";	        # R1 or R0

# from Invitrogen Superscript III and full-length cDNA (5' Cap) protocol
# my $adapter_5p="TCGTCGGGGACAACTTTGTACAAAAAAGTTGG";	# F1
# my $adapter_3p="CCCAACTTTCTTGTACAAAGTTGTCCCC";	# R1

### IMPORTANT !!! For cDNA PCR, adapter_5p is directly ligated onto the (3' end of the sense strand of)sequence (therefore the exact sequence is sequenced from 5' to 3'. BUT adapter_3p is used as a primer, so the actual sequence after sequencing is the reversed-complemented of adapter_3p!!!)

# 1. find "GTACTCTGCGTTGATACCACTGCTT" (rc___adapter_3p) and cut
# 2. determine strand :
#    2.1 yes upstream of "GTACTCTGCGTTGATACCACTGCTT" polyA (+), record polyA length, remove "AAGCAGTGGTATCAACGCAGAGTACATGGG"
#    2.2 no  upstream of "GTACTCTGCGTTGATACCACTGCTT" polyA (+), remove leading "AAGCAGTGGTATCAACGCAGAGTACTTTTT", record polyT length - 5, rc sequence
# 3. output substrings and polyA len (if length > minLne)


use strict;
my @time1 = localtime(time);
my @time2;

sub levenshtein($$){
    my @A=split //, lc shift;
    my @B=split //, lc shift;
    my @W=(0..@B);
    my ($i, $j, $cur, $next);
    for $i (0..$#A){
	$cur=$i+1;
	for $j (0..$#B){
	    $next=min( $W[$j+1]+1, $cur+1, ($A[$i] ne $B[$j])+$W[$j] );
            $W[$j]=$cur;
            $cur=$next;
	}
	$W[@B]=$next;
    }
    return $next;
}

sub countN($){
    my @A=split //, lc shift; 
    my $count=0;
    for(@A) {
	if($_ eq "n"){$count++;}
	elsif($_ eq "."){$count++;}
    }
    return $count;
}

sub min($$$){
    if ($_[0] < $_[2]){ pop @_; } else { shift @_; }
    return $_[0] < $_[1]? $_[0]:$_[1];
}

my $ADALEN3p=length($adapter_3p);	# R1
my $ADALEN5p=length($adapter_5p);	# F1
my $tmp=$adapter_3p;
$tmp=~tr/[ATCG]/[TAGC]/;
my $adapter_3p_rc=scalar reverse $tmp;
$tmp=$adapter_5p;
$tmp=~tr/[ATCG]/[TAGC]/;
my $adapter_5p_rc=scalar reverse $tmp;

my $kmer=6;
my $StepSize=3;
my %uniqF1rc;
my %uniqR1rc;
my %uniq_sanity;
for(my $i=0; $i<=($ADALEN5p-$kmer); $i=$i+$StepSize) {
    my $tmp_seq=substr($adapter_5p_rc,$i,$kmer);
    if (exists $uniq_sanity{$tmp_seq}) { $uniq_sanity{$tmp_seq}=$uniq_sanity{$tmp_seq}."\t".$i; }
    else { $uniq_sanity{$tmp_seq}=$i; }
}
foreach my $id (keys %uniq_sanity) {
    my @a=split("\t",$uniq_sanity{$id});
    if (scalar(@a) eq 1) {
	$uniqF1rc{$uniq_sanity{$id}}=$id;
    }
}

%uniq_sanity=();
for(my $i=0; $i<=($ADALEN3p-$kmer); $i=$i+$StepSize) {
    my $tmp_seq=substr($adapter_3p_rc,$i,$kmer);
    if (exists $uniq_sanity{$tmp_seq}) { $uniq_sanity{$tmp_seq}=$uniq_sanity{$tmp_seq}."\t".$i; }
    else { $uniq_sanity{$tmp_seq}=$i; }
}
foreach my $id (keys %uniq_sanity) {
    my @a=split("\t",$uniq_sanity{$id});
    if (scalar(@a) eq 1) {
	$uniqR1rc{$uniq_sanity{$id}}=$id;
    }
}

sub findpolyA ($) {
    my @A=split //, lc shift;
    my $Nr=scalar(@A) - 1;
    my $Cnt=0;
    my $window=10;
    my $len=$window;
    for (my $k=0; $k<=$len; $k++) { if ($A[$Nr - $k] eq "a") {$Cnt++;} }
    if ($Cnt/$len < 0.8) {return 0; }
    my $f=0;
    for (my $k=0; $k<=4; $k++) {  if ($A[$Nr - $len + $k] ne "a") {$f++;} else {last;} }
    while (($Cnt/$len >= 0.8) and ($f <= 3) and ($len < $Nr)) {
	$len+=$window/2;
	for (my $k=0; $k<=$len; $k++) { if ($A[$Nr - $k] eq "a") {$Cnt++;} }
	for (my $k=0; $k<=4; $k++) {  if ($A[$Nr - $len + $k] ne "a") {$f++;} else {last;} }
    }
    my $plen=$len - $f + 1;
    if ($plen < 5) { return 0; }
    else { return $plen; }
}

open IN1,$filein1;
if ($remakeID eq 1) {
	#open OUT1,">".$filein1.".t15";
	open OUT2,">".$filein1.".t15.info";
	#open OUT3,">".$filein1.".discard";
	#open OUTerr,">".$filein1.".t15_dist";
}
else {
	open OUT1,">".$filein1.".t15";
	#open OUT2,">".$filein1.".t15.info";
	open OUT3,">".$filein1.".discard";
	#open OUTerr,">".$filein1.".t15_dist";
}


my $pinfo="";
my $filetype;
my $total=0;
my $saved=0;
my $untrim=0;
my %leng;
while(<IN1>)
{
    my $id=$_;
    chomp $id;
    my $sequence=<IN1>;
    chomp $sequence;       

    my $reflen=length($sequence);
    my $quality="";
    my $quality2="";
    if($id=~m/^>/) {$filetype="fa";}
    elsif($id=~m/^@/) {$filetype="fq"; <IN1>; $quality=<IN1>; chomp $quality; }
    
    if ($reflen < $minLen) {
		if ($remakeID eq 0) { print OUT3 $id,"\n",$sequence,"\n"; }
		next;
    }
    
    my $processed=0;
    #if ($verbose) {print "3p_adapter\t",$adapter_3p,"\n","3p_adapter_rc\t",$adapter_3p_rc,"\n","5p_adapter\t",$adapter_5p,"\n","5p_adapter_rc\t",$adapter_5p_rc,"\n";}
    if ($verbose) {print  "\nStart of phase-1, looking for $adapter_3p_rc \n",$id,"\n",$sequence,"\t",$reflen,"\n";}
    $pinfo="";
    
    # 1. find "GTACTCTGCGTTGATACCACTGCTT" (rc___adapter_3p) and cut
    my $f=0;	
    my %hits;
    my $seq=$sequence;
    my %Subseq;
	my %CUT;
    my $CntSS=-1;
    my $SSstart=0;
    foreach my $anchor(sort {$a <=> $b} keys %uniqR1rc) {
        my $offset=-1;
        my $found=index($seq, $uniqR1rc{$anchor}, $offset);
        while($found != -1) {
            my $flag=0;
            for (my $p=1; ($p <= $f)and($flag == 0); $p++){
                my @info=split("\t",$hits{$p});
				if (abs(($found - $anchor) - ($info[1] - $info[0])) <= $kmer){ $info[2]++; $hits{$p}=join("\t",@info); $flag=1; }
            }					
            if ($flag eq 0) { $f++; $hits{$f}=join("\t",$anchor,$found,1); }
            $offset=$found+1;
            $found=index($seq, $uniqR1rc{$anchor}, $offset);
        }
    }
    if ($f ne 0) {
		my %result;
		for(my $i=1; $i<=$f; $i++ ) {
		    my @tmp=split("\t",$hits{$i});
		    if ($verbose) { print join("\t",@tmp),"\t",$i,"\t",$uniqR1rc{$tmp[0]},"\n"; }
		    $result{$tmp[1]}{$tmp[0]}{$tmp[2]}=join("\t",@tmp);
		}	
		foreach my $q_start (sort{$a <=> $b} keys %result) {
		    foreach my $start (sort{$a <=> $b} keys %{$result{$q_start}}) {
				foreach my $kmer (sort{$b <=> $a} keys %{$result{$q_start}{$start}}) {
				    my @tmp=split("\t",$result{$q_start}{$start}{$kmer});
				    my $offset=$tmp[1] > $tmp[0] ? ($tmp[1] - $tmp[0]) : 0;
				    my $query=substr($seq,$offset,$ADALEN3p);
				    my $qlen=length($query);
				    my $adapter_tmp=substr($adapter_3p_rc,0,$qlen);
				    my $dist=levenshtein($query,$adapter_tmp) - countN($query);
				    if ($verbose) { print join("\t",@tmp),"\n",$query,"\n",$adapter_tmp,"\n","edit-distance = ",$dist,"\n","offset = \t",$offset,"\n"; }
				    if ($dist/$qlen < $ErrorRate ) {
						if ($SSstart >= $reflen) {next;}
						my $tmp_seq=substr($sequence,$SSstart,($offset - $SSstart));
						$CntSS++;
						$Subseq{$CntSS}=$tmp_seq;
						$CUT{$CntSS}=1;
						$SSstart+=$offset;
						$SSstart+=$ADALEN3p;
						if ($verbose) { print "save : $tmp_seq \n"; }
					}
		        }
		    }
		}
    }
	if ($SSstart < $reflen) {
		$CntSS++;
		$CUT{$CntSS}=0;
		$Subseq{$CntSS}=substr($sequence,$SSstart);
	}
    undef %hits;
    # if ($CntSS eq -1) { $CntSS++; $Subseq{$CntSS}=$sequence; }
    if ($verbose) {
		print "============= status ============= \n";
		print  "phase-1 : end\n",$id,"\n",$sequence,"\t",$reflen,"\n";
		for(my $i=0; $i<=$CntSS; $i++) {
		    print ">".$i,"\n",$Subseq{$i},"\n";
		}
		print "============= status ============= \n";
    }    
    
    if ($verbose) {print  "\nStart of phase-2\n";}
    # 2. determine strand :
    #    2.1 yes upstream of "GTACTCTGCGTTGATACCACTGCTT" polyA (+), record polyA length, remove "AAGCAGTGGTATCAACGCAGAGTACATGGG"
    #    2.2 no  upstream of "GTACTCTGCGTTGATACCACTGCTT" polyA (+), remove leading "AAGCAGTGGTATCAACGCAGAGTACTTTTT", record polyT length - 5, rc sequence
    for(my $i=0; $i<=$CntSS; $i++) {
		my $seq=$Subseq{$i};
		my $reflen=length($seq);
		if ($reflen < $minLen) {
			if ($remakeID eq 0) {
				print OUT3 $id."___".$i,"\n",$seq,"\n";
			}
		    next;
		}
		# search from the 3' end for polyA
		my $polyAlen=findpolyA($seq);
		if ($polyAlen > 0) {
		    #    2.1 yes upstream of "GTACTCTGCGTTGATACCACTGCTT" polyA (+), record polyA length, remove "AAGCAGTGGTATCAACGCAGAGTACATGGG"
		    my $new_seq=substr($seq,0,($reflen - $polyAlen));
		    # find possible leading "AAGCAGTGGTATCAACGCAGAGTACATGGG"
		    my $mindiff=999;
		    my $minpos=-1;
		    for(my $j=0; $j<=100; $j++) {
				my $adapter_tmp=$adapter_5p;
				my $query=substr($new_seq,$j,$ADALEN5p);
				if (length($query) < $ADALEN5p) {last;}
				my $dist=levenshtein($query,$adapter_tmp) - countN($query);
				if ($dist/$ADALEN5p < $ErrorRate) {
				    if ($mindiff > $dist) {
						$mindiff=$dist;
						$minpos=$j;
				    }
				}
		    }
		    if ($minpos eq -1) {
				# B3
				if (length($new_seq) > $minLen) {
					if ($remakeID eq 0) { print OUT1 $id."_SUB".$i."_uNAp_".$polyAlen,"\n",$new_seq,"\n"; }
					else { print OUT2 $id."_SUB".$i."_".$polyAlen,"\t",$id."_SUB".$i."_uNAp_".$polyAlen,"\n"; }
				}
				if ($verbose) { print $id."_SUB".$i."_uNAp_".$polyAlen,"\n",$new_seq,"\n\n"; };
		    }
		    else {
				# B1
				my $nnew_seq=substr($new_seq,($minpos+$ADALEN5p));
				if (length($nnew_seq) > $minLen) {
					if ($remakeID eq 0) { print OUT1 $id."_SUB".$i."_pNAp_".$polyAlen,"\n",$nnew_seq,"\n"; }
					else { print OUT2 $id."_SUB".$i."_".$polyAlen,"\t",$id."_SUB".$i."_pNAp_".$polyAlen,"\n"; }
				}
				if ($verbose) { print $id."_SUB".$i."_pNAp_".$polyAlen,"\n",$nnew_seq,"\n\n"; };
		    }
		}
		else {
		    #    2.2 no  upstream of "GTACTCTGCGTTGATACCACTGCTT" polyA (+), remove leading "AAGCAGTGGTATCAACGCAGAGTACTTTTT", record polyT length - 5, rc sequence
		    my $tmp=$seq;
		    $tmp=~tr/[ATCG]/[TAGC]/;
		    my $new_seq=scalar reverse $tmp;
		    $reflen=length($new_seq);
		    if ($verbose) { print $new_seq,"\t",$reflen,"\n"; }
		    # find possible ending "GTACTCTGCGTTGATACCACTGCTT"
		    my $mindiff=999;
		    my $minpos=-1;
		    for(my $j=0; $j<=100; $j++) {
				my $adapter_tmp=$adapter_3p_rc;
				my $query=substr($new_seq,(0 - $j - $ADALEN3p),$ADALEN3p);
				if (length($query) < $ADALEN3p) {last;}
				if ($verbose) { print "adapt : ",$adapter_tmp,"\n","query : ",$query,"\n"; }
				my $dist=levenshtein($query,$adapter_tmp) - countN($query);
				if ($dist/$ADALEN3p < $ErrorRate) {
				    if ($mindiff > $dist) {
					$mindiff=$dist;
					$minpos=$j;
				    }
				}
		    }
		    if ($verbose) { print "mindiff = ",$mindiff,"\tminpos = ", $minpos,"\n"; }
		    if ($minpos eq -1) {
				# A3
				$polyAlen=findpolyA($new_seq);
				if ($polyAlen > 0) {
					my $nnew_seq=substr($new_seq,0,($reflen - $polyAlen));
					if ($CUT{$i} eq 1) {
						if (length($nnew_seq) > $minLen) {
							if ($remakeID eq 0) { print OUT1 $id."_SUB".$i."_uTNp_".$polyAlen,"\n",$nnew_seq,"\n"; }
							else { print OUT2 $id."_SUB".$i."_".$polyAlen,"\t",$id."_SUB".$i."_uTNp_".$polyAlen,"\n" }
						}
						if ($verbose) { print $id."_SUB".$i."_uTNp_".$polyAlen,"\n",$nnew_seq,"\n"; };
					}
					else {
						if (length($nnew_seq) > $minLen) {
							if ($remakeID eq 0) { print OUT1 $id."_SUB".$i."_uTNu_".$polyAlen,"\n",$nnew_seq,"\n"; }
							else {print OUT2 $id."_SUB".$i."_".$polyAlen,"\t",$id."_SUB".$i."_uTNu_".$polyAlen,"\n";}
						}
						if ($verbose) { print $id."_SUB".$i."_uTNu_".$polyAlen,"\n",$nnew_seq,"\n"; };
					}
				}
				else {
					if ($CUT{$i} eq 1) {
						if (length($new_seq) > $minLen) {
							if ($remakeID eq 0) { print OUT1 $id."_SUB".$i."_utNp_".$polyAlen,"\n",$new_seq,"\n"; }
							else { print OUT2 $id."_SUB".$i."_".$polyAlen,"\t",$id."_SUB".$i."_utNp_".$polyAlen,"\n";}
						}
						if ($verbose) { print $id."_SUB".$i."_utNp_".$polyAlen,"\n",$new_seq,"\n"; };
					}
					else {
						if (length($new_seq) > $minLen) {
							if ($remakeID eq 0) { print OUT1 $id."_SUB".$i."_utNu_".$polyAlen,"\n",$new_seq,"\n"; }
							else { print OUT2 $id."_SUB".$i."_".$polyAlen,"\t",$id."_SUB".$i."_utNu_".$polyAlen,"\n"; }
						}
						if ($verbose) { print $id."_SUB".$i."_utNu_".$polyAlen,"\n",$new_seq,"\n"; };
					}
				}
		    }
		    else {
				my $nnew_seq=substr($new_seq,0,($reflen - $minpos - $ADALEN5p));
				$reflen=length($nnew_seq);
				$polyAlen=findpolyA($nnew_seq);
				if ($polyAlen > 0) {
					# A1 + A2
				    my $nnnew_seq=substr($nnew_seq,0,($reflen - $polyAlen));
					if ($CUT{$i} eq 1) {
						if (length($nnnew_seq) > $minLen) {
							if ($remakeID eq 0) { print OUT1 $id."_SUB".$i."_pTNp_".$polyAlen,"\n",$nnnew_seq,"\n"; }
							else { print OUT2 $id."_SUB".$i."_".$polyAlen,"\t",$id."_SUB".$i."_pTNp_".$polyAlen,"\n"; }
						}
						if ($verbose) { print $id."_SUB".$i."_pTNp_".$polyAlen,"\n",$nnnew_seq,"\n"; };
					}
					else {
						if (length($nnnew_seq) > $minLen) {
							if ($remakeID eq 0) { print OUT1 $id."_SUB".$i."_pTNu_".$polyAlen,"\n",$nnnew_seq,"\n"; }
							else { print OUT2 $id."_SUB".$i."_".$polyAlen,"\t",$id."_SUB".$i."_pTNu_".$polyAlen,"\n"; }
						}
						if ($verbose) { print $id."_SUB".$i."_pTNu_".$polyAlen,"\n",$nnnew_seq,"\n"; };
					}
				}
				else {
					# B2
					my $tmp_seq=$nnew_seq;
					$tmp_seq=~tr/[ATCG]/[TAGC]/;
					my $ttmp_seq=scalar reverse $tmp_seq;
					if ($CUT{$i} eq 1) {
						if (length($ttmp_seq) > $minLen) {
							if ($remakeID eq 0) {  print OUT1 $id."_SUB".$i."_pNAu_".$polyAlen,"\n",$ttmp_seq,"\n"; }
							else { print OUT2 $id."_SUB".$i."_".$polyAlen,"\t",$id."_SUB".$i."_pNAu_".$polyAlen,"\n"; }
						}
						if ($verbose) { print $id."_SUB".$i."_pNAu_".$polyAlen,"\n",$ttmp_seq,"\n"; };
					}
					else {
						if (length($ttmp_seq) > $minLen) {
							if ($remakeID eq 0) {  print OUT1 $id."_SUB".$i."_pNau_".$polyAlen,"\n",$ttmp_seq,"\n"; }
							else { print OUT2 $id."_SUB".$i."_".$polyAlen,"\t",$id."_SUB".$i."_pNau_".$polyAlen,"\n"; }
						}
						if ($verbose) { print $id."_SUB".$i."_pNau_".$polyAlen,"\n",$ttmp_seq,"\n"; };
					}
				    
				}
		    }
		}
    }
    # 3. output substrings and polyA len (if length > minLne)
}
#							CUT		polyAlen	cut2	polyTlen	Code
# A1 : -> TTTNNNNN <-		1		0			1		1			pTNp
# A2 : -> TTTNNN   			0		0			1		1			pTNu

# A3 :      TNNNNN <-		1		0			0		0			uTNp / uTNu / utNp / utNu
# B2 : -> NNNNNA  			0		0			1		0			pNau

# B1 : -> NNNNNAAA <-		1		1			1					pNAp
# B3 :      NNNAAA <-		1		1			0					uNAp

