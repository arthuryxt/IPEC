#!/usr/bin/perl -w
use strict;
die "Usage: $0   \"input_refFlat\"    \"input_wiggle_list\"    \"output_basename\"   \"window_size==100\"    \"cutoff==0.05\"   \"mincov==5\"   \(strandness==-1\)" if (@ARGV < 3);
my $filein1=$ARGV[0];
my $filein2=$ARGV[1];   # note that assume all reads are on the opposite strand of the transcripts! # make sure the wiggle file contains ONLY ONE "minus" or "plus" in the track_def line
my $fileout=$ARGV[2];
my $window=100;
if (scalar(@ARGV) > 3) { $window=$ARGV[3]; }
my $cutoff=0.05;
if (scalar(@ARGV) > 4) { $cutoff=$ARGV[4]; }
my $mincov=1;
if (scalar(@ARGV) > 5) { $mincov=$ARGV[5]; }
my $strandness=0;    # sequensing_opposite_strand_of_transcript==-1; non-strandess==0; sequensing_same_strand_of_transcript==1;
if (scalar(@ARGV) > 6) { $strandness=$ARGV[6]; }
if (($strandness ne -1) and ($strandness ne 0) and ($strandness ne 1)) { die "strandness must be one of [-1, 0, 1] \n ";}
my $step=int($window/2);


my %uniq;
my %lastExon;
open(IN1, $filein1) or die "Cannot open input_refFlat file : $filein1";
while (<IN1>) {
    chomp;
    next if (m/^@/);
    next if (m/^#/);
    next if (m/^track/);
    my @a=split("\t",$_);
    if ($a[2]!~m/^chr/) { $a[2]="chr".$a[2]; }
    #     chr    strand GL    transcript
    $uniq{$a[2]}{$a[3]}{$a[0]}{$a[1]}=join("\t",$_);
    if ($a[3] eq "+") {
        my @Start=split(/\,/,$a[9]);
        my @End=split(/\,/,$a[10]);
        if (exists $lastExon{$a[2]}{$a[3]}{$a[0]}) {
            my @tmp=split("\t",$lastExon{$a[2]}{$a[3]}{$a[0]});
            if ($tmp[1] >= $End[-1]) {
                # the recorded was longer
                if ($tmp[0] < $Start[-1]) { $lastExon{$a[2]}{$a[3]}{$a[0]}=$Start[-1]."\t".$tmp[1]; }
            }
            else {
                if ($tmp[0] > $Start[-1]) { $lastExon{$a[2]}{$a[3]}{$a[0]}=$tmp[0]."\t".$End[-1]; }
                else { $lastExon{$a[2]}{$a[3]}{$a[0]}=$Start[-1]."\t".$End[-1]; }
            }
        }
        else {
            $lastExon{$a[2]}{$a[3]}{$a[0]}=$Start[-1]."\t".$End[-1];
        }
    }
    else {
        my @Start=split(/\,/,$a[9]);
        my @End=split(/\,/,$a[10]);
        if (exists $lastExon{$a[2]}{$a[3]}{$a[0]}) {
            my @tmp=split("\t",$lastExon{$a[2]}{$a[3]}{$a[0]});
            if ($tmp[0] <= $Start[0]) {
                # the recorded was longer
                if ($tmp[1] > $End[0]) { $lastExon{$a[2]}{$a[3]}{$a[0]}=$tmp[0]."\t".$End[0]; }
            }
            else {
                if ($tmp[1] < $End[0]) { $lastExon{$a[2]}{$a[3]}{$a[0]}=$Start[0]."\t".$tmp[1]; }
                else { $lastExon{$a[2]}{$a[3]}{$a[0]}=$Start[0]."\t".$End[0]; }
            }
        }
        else {
            $lastExon{$a[2]}{$a[3]}{$a[0]}=$Start[0]."\t".$End[0];
        }
    }
}
close IN1;

my %checked;
open(OUT1, ">".$fileout.".FLT");
open(OUT2, ">".$fileout.".nonflt");
open(OUT3, ">".$fileout.".3endinfo");
open(IN2, $filein2) or die "Cannot open wiggle file : $filein2";
while (<IN2>) {
    chomp;
    my $file=$_;
    open(INw, $file) or die "Cannot open wiggle file : $file";
    print "reading wiggle file : $file \n";
    my $strand="+";
    my $chr="";
    my %cov;
    while (<INw>) {
        chomp;
        if (m/^track/) {
            if ($_=~m/plus/) {
                if ($strandness eq -1) { $strand="-"; }
                elsif ($strandness eq 1) { $strand="+"; }
                else { $strand="-"; }
            }
            elsif ($_=~m/minus/) {
                if ($strandness eq -1) { $strand="+"; }
                elsif ($strandness eq 1) { $strand="-"; }
                else { $strand="+"; }
            }
            else {
                print "wiggle file without strand-information : $file \n";
                last;
            }
        }
        else {
            my @a=split("\t",$_);
            if ($chr eq "") {
                $chr=$a[0];
                for(my $i=$a[1]; $i<$a[2]; $i++) { $cov{$chr}{$i}=int(1+$a[3]); }
                print "reading $chr $strand \n";
            }
            elsif ($chr eq $a[0]) {
                for(my $i=$a[1]; $i<$a[2]; $i++) { $cov{$chr}{$i}=int(1+$a[3]); }
            }
            else {
                # encouter a new chromosome, report, delete hash, and restart
                print "current chr = $chr, new chr = $a[0],\n";
                if ($strand eq "+") {
                    foreach my $loci (keys %{$lastExon{$chr}{$strand}}) {
                        my @exon=split("\t",$lastExon{$chr}{$strand}{$loci});
                        my $CL=0;
                        my $len=$exon[1] - $exon[0];
                        for(my $i=$exon[0]; $i<=$exon[1]; $i++) {
                            if (exists $cov{$chr}{$i}) { $CL+=$cov{$chr}{$i}; }
                        }
                        $CL=sprintf("%.3f",$CL/$len);
                        my @covarray;
                        $covarray[0]=$CL;
                        my $ext=0;
                        while ($covarray[$ext] > $cutoff * $CL) {
                            $ext++;
                            my $tmpcov=0;
                            for(my $i=0; $i<$window; $i++) {
                                if (exists $cov{$chr}{($exon[1] - 50 + 50*$ext + $i)}) { $tmpcov+=$cov{$chr}{($exon[1] - 50 + 50*$ext + $i)}; }
                            }
                            $covarray[$ext]=sprintf("%.3f",$tmpcov/$window);
                        }
                        print OUT3 join("\t",$chr,$exon[0],$exon[1],$loci,$ext,$strand,join(",",@covarray)),"\n";
                        my $drop=0;
                        for(my $i=0; $i<=$ext; $i++) { if($covarray[$i] < $mincov){ $drop=$i; last; } }
                        if ($drop <= 2) {
                            # true FLT
                            foreach my $tid (keys %{$uniq{$chr}{$strand}{$loci}}) {
                                print OUT1 $uniq{$chr}{$strand}{$loci}{$tid},"\n";
                            }
                        }
                        else {
                            foreach my $tid (keys %{$uniq{$chr}{$strand}{$loci}}) {
                                print OUT2 $uniq{$chr}{$strand}{$loci}{$tid},"\n";
                            }
                        }
                    }
                }
                else {
                    foreach my $loci (keys %{$lastExon{$chr}{$strand}}) {
                        my @exon=split("\t",$lastExon{$chr}{$strand}{$loci});
                        my $CL=0;
                        my $len=$exon[1] - $exon[0];
                        for(my $i=$exon[0]; $i<=$exon[1]; $i++) {
                            if (exists $cov{$chr}{$i}) { $CL+=$cov{$chr}{$i}; }
                        }
                        $CL=sprintf("%.3f",$CL/$len);
                        my @covarray;
                        $covarray[0]=$CL;
                        my $ext=0;
                        while ($covarray[$ext] > $cutoff * $CL) {
                            $ext++;
                            my $tmpcov=0;
                            for(my $i=0; $i<$window; $i++) {
                                if (exists $cov{$chr}{($exon[0] + 50 - 50*$ext - $i)}) { $tmpcov+=$cov{$chr}{($exon[0] + 50 - 50*$ext - $i)}; }
                            }
                            $covarray[$ext]=sprintf("%.3f",$tmpcov/$window);
                        }
                        print OUT3 join("\t",$chr,$exon[0],$exon[1],$loci,$ext,$strand,join(",",@covarray)),"\n";
                        my $drop=0;
                        for(my $i=0; $i<=$ext; $i++) { if($covarray[$i] < $mincov){ $drop=$i; last; } }
                        if ($drop <= 2) {
                            # true FLT
                            foreach my $tid (keys %{$uniq{$chr}{$strand}{$loci}}) {
                                print OUT1 $uniq{$chr}{$strand}{$loci}{$tid},"\n";
                            }
                        }
                        else {
                            foreach my $tid (keys %{$uniq{$chr}{$strand}{$loci}}) {
                                print OUT2 $uniq{$chr}{$strand}{$loci}{$tid},"\n";
                            }
                        }
                    }
                }
                
                
                # delete hash
                %cov=();
                # restart
                $chr=$a[0];
                for(my $i=$a[1]; $i<$a[2]; $i++) { $cov{$chr}{$i}=int(1+$a[3]); }
                print "reading $chr $strand \n";
            }
        }
    }
    {
            print "processing $chr $strand\n";
                if ($strand eq "+") {
                    foreach my $loci (keys %{$lastExon{$chr}{$strand}}) {
                        my @exon=split("\t",$lastExon{$chr}{$strand}{$loci});
                        my $CL=0;
                        my $len=$exon[1] - $exon[0];
                        for(my $i=$exon[0]; $i<=$exon[1]; $i++) {
                            if (exists $cov{$chr}{$i}) { $CL+=$cov{$chr}{$i}; }
                        }
                        $CL=sprintf("%.3f",$CL/$len);
                        my @covarray;
                        $covarray[0]=$CL;
                        my $ext=0;
                        while ($covarray[$ext] > $cutoff * $CL) {
                            $ext++;
                            my $tmpcov=0;
                            for(my $i=0; $i<$window; $i++) {
                                if (exists $cov{$chr}{($exon[1] - 50 + 50*$ext + $i)}) { $tmpcov+=$cov{$chr}{($exon[1] - 50 + 50*$ext + $i)}; }
                            }
                            $covarray[$ext]=sprintf("%.3f",$tmpcov/$window);
                        }
                        print OUT3 join("\t",$chr,$exon[0],$exon[1],$loci,$ext,$strand,join(",",@covarray)),"\n";
                        my $drop=0;
                        for(my $i=0; $i<=$ext; $i++) { if($covarray[$i] < $mincov){ $drop=$i; last; } }
                        if ($drop <= 2) {
                            # true FLT
                            foreach my $tid (keys %{$uniq{$chr}{$strand}{$loci}}) {
                                print OUT1 $uniq{$chr}{$strand}{$loci}{$tid},"\n";
                            }
                        }
                        else {
                            foreach my $tid (keys %{$uniq{$chr}{$strand}{$loci}}) {
                                print OUT2 $uniq{$chr}{$strand}{$loci}{$tid},"\n";
                            }
                        }
                    }
                }
                else {
                    foreach my $loci (keys %{$lastExon{$chr}{$strand}}) {
                        my @exon=split("\t",$lastExon{$chr}{$strand}{$loci});
                        my $CL=0;
                        my $len=$exon[1] - $exon[0];
                        for(my $i=$exon[0]; $i<=$exon[1]; $i++) {
                            if (exists $cov{$chr}{$i}) { $CL+=$cov{$chr}{$i}; }
                        }
                        $CL=sprintf("%.3f",$CL/$len);
                        my @covarray;
                        $covarray[0]=$CL;
                        my $ext=0;
                        while ($covarray[$ext] > $cutoff * $CL) {
                            $ext++;
                            my $tmpcov=0;
                            for(my $i=0; $i<$window; $i++) {
                                if (exists $cov{$chr}{($exon[0] + 50 - 50*$ext - $i)}) { $tmpcov+=$cov{$chr}{($exon[0] + 50 - 50*$ext - $i)}; }
                            }
                            $covarray[$ext]=sprintf("%.3f",$tmpcov/$window);
                        }
                        print OUT3 join("\t",$chr,$exon[0],$exon[1],$loci,$ext,$strand,join(",",@covarray)),"\n";
                        my $drop=0;
                        for(my $i=0; $i<=$ext; $i++) { if($covarray[$i] < $mincov){ $drop=$i; last; } }
                        if ($drop <= 2) {
                            # true FLT
                            foreach my $tid (keys %{$uniq{$chr}{$strand}{$loci}}) {
                                print OUT1 $uniq{$chr}{$strand}{$loci}{$tid},"\n";
                            }
                        }
                        else {
                            foreach my $tid (keys %{$uniq{$chr}{$strand}{$loci}}) {
                                print OUT2 $uniq{$chr}{$strand}{$loci}{$tid},"\n";
                            }
                        }
                    }
                }
    }
    close INw;
    $strand="+";
    $chr="";
}


