#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"config_file\"   \"output_sh\"  " if (@ARGV < 2);
my $filein=$ARGV[0];
my $fileout=$ARGV[1];
my %SPEC;
open(IN,$filein) or die "Cannot open config_file $filein";
while (<IN>) {
    chomp;
    my @a=split("\t",$_);
    if ((scalar(@a) < 2) or ($a[0] eq "")) { next; }
    $SPEC{$a[0]}=$a[1];
}
close IN;
my $command="rm -rf $fileout";
system($command);

my $threads=1;
my $remove_temp="yes";
my $output_folder="output";
# check if all parameters are set
if (!exists $SPEC{"PacBio_fasta"})        { die "ERROR: A fasta file of PacBio read must by specified in the config_file $filein";}
if (!exists $SPEC{"Illumina_fasta_list"}) { die "ERROR: A list of Illumina reads must by specified in the config_file $filein";}
# the file needs to contain FULL path to the files.
# if the files are "NGSa","NGSb",...,"NGSz", one can prepare the file by : ls -d $PWD/NGS? > NGS_file_list.txt
if (!exists $SPEC{"Bowtie2_folder"})      { die "ERROR: The folder containing Bowtie2 executables must by specified in the config_file $filein";}
if (!exists $SPEC{"IPEC_folder"})         { die "ERROR: The folder containing ipec scripts must by specified in the config_file $filein";}
if (!exists $SPEC{"Bowtie2_folder"})      { $SPEC{"PacBio_fasta"}=$SPEC{"PacBio_fasta"}}
if (!exists $SPEC{"Illumina_fasta_list"}) { $SPEC{"Illumina_fasta_list"}=$SPEC{"Illumina_fasta_list"}}
if (!exists $SPEC{"Bowtie2_folder"})      { $SPEC{"Bowtie2_folder"}=$SPEC{"Bowtie2_folder"}}
if (!exists $SPEC{"IPEC_folder"})         { $SPEC{"IPEC_folder"}=$SPEC{"IPEC_folder"}}
if (exists $SPEC{"Bowtie2_threads"})      { $threads=$SPEC{"Bowtie2_threads"}; }
if (exists $SPEC{"Remove_tmp_files"})     { $remove_temp=$SPEC{"Remove_tmp_files"}; }
if (exists $SPEC{"Output_folder"})        { $output_folder=$SPEC{"Output_folder"}; }


open(OUT, ">".$fileout) or die "ERROR: Cannot open output_sh file $fileout";
print OUT "#!/bin/bash\n\n";
print OUT "date\n";
print OUT "#Step1\n";
print OUT "echo \"Step1 : Build index for PacBio reads. Starting...\" \n";
$command="mkdir $output_folder";
print OUT $command,"\n";
$command="cd $output_folder";
print OUT $command,"\n";
$command="mkdir \"1_Index\"";
print OUT $command,"\n";
$command="mkdir \"2_Mapping\"";
print OUT $command,"\n";
$command="mkdir \"3_EC\"";
print OUT $command,"\n";
$command="mkdir \"Results\"";
print OUT $command,"\n";
$command=$SPEC{"Bowtie2_folder"}."/bowtie2-build -o 1 --threads ".$threads." ".$SPEC{"PacBio_fasta"}." ./1_Index/PacBio";
print OUT $command,"\n";
print OUT "echo \"Step1 : Build index for PacBio reads. Finished\" \n\n\n";


print OUT "#Step2\n";
print OUT "date\n";
print OUT "echo \"Step2 : Mapping Illumina reads to PacBio reads using sensitive local alignment. Starting...\" \n";
open IN,$SPEC{"Illumina_fasta_list"};
my $filecnt=0;
while(<IN>){
    chomp;
    my $input=$_;
    my $output="NGS_".$filecnt;
    $filecnt++;
    $command=$SPEC{"Bowtie2_folder"}."bowtie2 -p ".$threads." -f --phred64 --local -D 40, -R 3, -N 1 -L 20 -i S,1,0.50 -k 5000 --no-head --ma 2 --mp 6 --rdg 3,3 --rfg 3,3 -x ./1_index/PacBio ".$input." > ./2_Mapping/".$output.".bt2aln";
    print OUT $command,"\n";
}
close IN;
print OUT "echo \"Step2 : Mapping Illumina reads to PacBio reads using sensitive local alignment. Finished\" \n\n\n";


print OUT "#Step3\n";
print OUT "date\n";
print OUT "echo \"Step3 : Processing  correction using crowd wisdom. Starting...\" \n";
for(my $i=0; $i<$filecnt; $i++) {
    my $input="NGS_".$i.".bt2aln";
    my $output="NGS_".$i.".fa";
    $command="perl ".$SPEC{"IPEC_folder"}."/IPEC_1.pl ".$SPEC{"PacBio_fasta"}." ./2_Mapping/".$input." ./3_EC/".$output." 2";
    print OUT $command,"\n";
}
print OUT "echo \"Step3 : Processing correction using crowd wisdom. Finished\" \n\n\n";


print OUT "#Step4\n";
print OUT "date\n";
print OUT "echo \"Step4 : Merge information from partially (NGS-)corrected into fully corrected. Starting...\" \n";
$command="ls ./3_EC/NGS_\*fa.mod \| sort -u \> ipec.mod.filelist";
print OUT $command,"\n";
$command="perl ".$SPEC{"IPEC_folder"}."/IPEC_1_merge.pl ./Results/corrected ipec.mod.filelist";
print OUT $command,"\n";
print OUT "echo \"Step4 : Merge information from partially (NGS-)corrected into fully corrected. Finished\" \n\n\n";

if($remove_temp eq "yes") {
    $command="rm -rf ./1_Index";
    print OUT $command,"\n";
    $command="rm -rf ./2_Mapping";
    print OUT $command,"\n";
    $command="rm -rf ./3_EC";
    print OUT $command,"\n";
    $command="rm -rf ipec.mod.filelist";
    print OUT $command,"\n";
}

print OUT "date\n";

close OUT;

