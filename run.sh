#!/bin/bash

date
#Step1
echo "Step1 : Build index for PacBio reads. Starting..." 
mkdir output
cd output
mkdir "1_Index"
mkdir "2_Mapping"
mkdir "3_EC"
mkdir "Results"
/Users/arthur/Desktop/bowtie2//bowtie2-build -o 1 --threads 2 /Users/arthur/Desktop/IPEC/data/PacBio.fa ./1_Index/PacBio
echo "Step1 : Build index for PacBio reads. Finished" 


#Step2
date
echo "Step2 : Mapping Illumina reads to PacBio reads using sensitive local alignment. Starting..." 
/Users/arthur/Desktop/bowtie2/bowtie2 -p 2 -f --phred64 --local -D 40, -R 3, -N 1 -L 20 -i S,1,0.50 -k 5000 --no-head --ma 2 --mp 6 --rdg 3,3 --rfg 3,3 -x ./1_index/PacBio /Users/arthur/Desktop/IPEC/data/NGSa > ./2_Mapping/NGS_0.bt2aln
/Users/arthur/Desktop/bowtie2/bowtie2 -p 2 -f --phred64 --local -D 40, -R 3, -N 1 -L 20 -i S,1,0.50 -k 5000 --no-head --ma 2 --mp 6 --rdg 3,3 --rfg 3,3 -x ./1_index/PacBio /Users/arthur/Desktop/IPEC/data/NGSb > ./2_Mapping/NGS_1.bt2aln
echo "Step2 : Mapping Illumina reads to PacBio reads using sensitive local alignment. Finished" 


#Step3
date
echo "Step3 : Processing  correction using crowd wisdom. Starting..." 
perl /Users/arthur/Desktop/IPEC/scripts//IPEC_1.pl /Users/arthur/Desktop/IPEC/data/PacBio.fa ./2_Mapping/NGS_0.bt2aln ./3_EC/NGS_0.fa 2
perl /Users/arthur/Desktop/IPEC/scripts//IPEC_1.pl /Users/arthur/Desktop/IPEC/data/PacBio.fa ./2_Mapping/NGS_1.bt2aln ./3_EC/NGS_1.fa 2
echo "Step3 : Processing correction using crowd wisdom. Finished" 


#Step4
date
echo "Step4 : Merge information from partially (NGS-)corrected into fully corrected. Starting..." 
ls ./3_EC/NGS_*fa.mod | sort -u > ipec.mod.filelist
perl /Users/arthur/Desktop/IPEC/scripts//IPEC_1_merge.pl ./Results/corrected ipec.mod.filelist
echo "Step4 : Merge information from partially (NGS-)corrected into fully corrected. Finished" 


rm -rf ./1_Index
rm -rf ./2_Mapping
rm -rf ./3_EC
rm -rf ipec.mod.filelist
date
