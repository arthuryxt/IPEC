#_1_# Prepare Fasta file for PacBio reads
It should be a FASTA file, and we refer to "PacBio.fa" in this tutorial.

#_2_# Prepare Fasta files for Illumina reads and create a file list
As the throughput of NGS grows, using ALL reads in one lane for correction can lead to memory crisis. We can split one big Illumina read file into several smaller ones, perform the correction with the smaller file, and later merge the voting information for a final version.
Let say, we have split the big Illumina file into 26 smaller ones: <NGSa>, <NGSb>, ..., <NGSz>.
We can create a file list with the following command:
ls -d $PWD/NGS? > NGS_file_list.txt

#_3_# Prepare parameters
The <SPEC.txt> file is a tab-delimited file, with the first column for keywords recognized by ipec tools, and the second column for user defined files/parameters.
We recommend to generate this file with the following perl script:
perl -e 'print join("\t","PacBio_fasta","/Users/arthur/Desktop/IPEC/data/PacBio.fa"),"\n";' > SPEC.txt
perl -e 'print join("\t","Illumina_fasta_list","/Users/arthur/Desktop/IPEC/data/NGS_file_list.txt"),"\n";' >> SPEC.txt
perl -e 'print join("\t","Bowtie2_folder","/Users/arthur/Desktop/bowtie2/"),"\n";' >> SPEC.txt
perl -e 'print join("\t","IPEC_folder","/Users/arthur/Desktop/IPEC/scripts/"),"\n";' >> SPEC.txt
perl -e 'print join("\t","Output_folder","output"),"\n";' >> SPEC.txt
perl -e 'print join("\t","Bowtie2_threads","2"),"\n";' >> SPEC.txt
perl -e 'print join("\t","Remove_tmp_files","yes"),"\n";' >> SPEC.txt

#_4_# Prepare pipeline bash file
perl /Users/arthur/Desktop/IPEC/scripts/ipec_MAKE.pl SPEC.txt run.sh

#_5_# Run pipeline
bash run.sh

#_RESULT_#
<output/Results/corrected.fa> is the corrected fasta
Comparing the raw and post-corrected sequences using BLAT to rn6.0:
browser details PacBio.fa     2204     5  3111  3113    91.0%  chr2   +    18392251  18416979  24729
browser details corrected.fa  2827     5  2840  2842   100.0%  chr2   +    18392251  18416979  24729
If you think another iteration is needed, simply use this file as "PacBio_fasta", re-make the <SPEC.txt> file, re-make the "run.sh" and run ipec again.










#_other_scripts_#
always keep the <reads_of_insert.fasta> files from each SMRT analysis run
to extract fasta files (one for ccs.fa, the other for subread.fa) from SMRT runs with two movies (old, pre-2013), use this: Pacbio_get_fasta.pl
to extract fasta files (one for ccs.fa, the other for subread.fa) from SMRT runs with three movies (post-2013),   use this: Pacbio2_get_fasta.pl
to remove cDNA primers (SMARTer protocol) from the subreads:
    perl Pacbio_trim_cDNA_primer2.pl subread.fa AAGCAGTGGTATCAACGCAGAGTAC AAGCAGTGGTATCAACGCAGAGTACATGGG 0.25 300
to seperate wells with CCS reads with those without:
    perl Pacbio_seperate_CCS_2.1.pl subread.fa.t15 reads_of_insert.fasta SubRead.fa.t15 CCS.fa mySample_1
to remove cDNA primers (SMARTer protocol) from the ccs:
    perl Pacbio_trim_cDNA_primer2.pl CCS.fa     AAGCAGTGGTATCAACGCAGAGTAC AAGCAGTGGTATCAACGCAGAGTACATGGG 0.25 300
    <SubReads.fa.t15> and <CCS.fa.t15> are the files with trimmed reads, with >= 300 bases, and the adapters were recognized with at most 25% error rate.
to map the pacbio reads to genome using GMAP:
    gmap_build -D /genome/gmapdb/ -d rn6 -k 12 -q 1 -s numeric-alpha /genome/rn6/rn6.fa
    gmap -D /genome/gmapdb/ -t 20 -d rn6 -k 12 --basesize 12 --sampling 1 -B 2 --min-intronlength 20 --canonical-mode 2 -f gff3_gene -n 0 -O --fails-as-input CCS.fa.t15 > CCS.gtf
to convert the gtf format into refFlat format:
    perl gtf2refFlat_gmap.pl CCS.gtf CCS.refFlat 1 /genome/gmapdb/rat/idlist
to select the representative transcript isoform per gene loci via "clustering" or "collapsing" (full-length) pacbio reads, with known exon-junctions from refseq and ensembl annotation, no supporting CAGE data, requiring at least 2 reads from illumina to support the exon-junctions
    perl cluster_refFlat.pl rn5_refGene_20140102.refFlat rn5_refGene_20140102.reFFlat no no no no 2
    perl cluster_refFlat.pl rn5_Ensembl_20140409.refFlat rn5_Ensembl_20140409.reFFlat no no no no 2
    cat rn5_refGene_20140102.refFlat | awk '$9 > 1' | perl -ne 'chomp; my @a=split("\t",$_); my @S=split(/\,/,$a[9]); my @E=split(/\,/,$a[10]); print join("\t",$a[2],$a[4],$a[6],"junc","1",$a[3],$a[4],$a[6],join(",",255,0,0),2,join(",",10,10),join(",",0,($a[6]-$a[4]-10))),"\n";' > rn5_refGene_20140102_junc.bed
    cat rn5_Ensembl_20140409.refFlat | awk '$9 > 1' | perl -ne 'chomp; my @a=split("\t",$_); my @S=split(/\,/,$a[9]); my @E=split(/\,/,$a[10]); print join("\t",$a[2],$a[4],$a[6],"junc","1",$a[3],$a[4],$a[6],join(",",255,0,0),2,join(",",10,10),join(",",0,($a[6]-$a[4]-10))),"\n";' > rn5_Ensembl_20140409_junc.bed
    cat rn5_refGene_20140102_junc.bed rn5_Ensembl_20140409_junc.bed  > rat201404_Ref_Ens_junction_bed
    perl cluster_refFlat.pl rat_pacbio.refFlat rat_pacbio_final rat_all_refFlat_cov rat_all_strand rat201404_Ref_Ens_junction_bed no 2

