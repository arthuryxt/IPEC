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
