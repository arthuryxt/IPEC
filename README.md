# IPEC

Tookit for Illumina-assisted PacBio Error Correction (IPEC).

* This is part of the rat hippocampus transcriptome annotation project, where hybrid sequencing technologies (PacBio + Illumina) were applied to get full-length transcriptome landscape for rat hippocampus.  


## Content

  1. PacBio-reads Error Correction


## Installation
1. Operating system requirement: any OS running perl (>= 5)
2. Download the toolkit: 
```git clone https://github.com/arthuryxt/IPEC.git```
3. Changed to the IPEC directory and run the scripts: 
```cd IPEC```
4. Download Bowtie2 executables if you don't yet have them. Supported version 2.3.5.1 (tested) and lower versions.

Expected time: a few seconds

## Usage

* __ipec_MAKE.pl__   
  _Description_: Make a bash pipeline to use Illumia reads for correcting sequencing errors in PacBio reads. The pipeline consists of 4 steps:
  1. Build bowtie2 index from fasta sequence of PacBio reads. As correction is done for each PacBio reads independently, one can limit the number of PacBio reads for one batch and make full use of parallelization.
  2. Mapping Illumina reads to PacBio reads using sensitive local alignment using Bowtie2.
  3. Correctiong PacBio reads using crowd wisdom (voting).
  4. Merge information from partially (NGS-)corrected into fully corrected reads, and report additional stats files.
  
  _Usage_: perl ipec_MAKE.pl <SPEC.txt> <run.sh> 



## Examples  

* __ipec_MAKE.pl__  
```
perl scripts/ipec_MAKE.pl SPEC.txt run.sh   
```

* __run.sh__   
```
bash run.sh 
```

Expected time: a couple of seconds for the toy example








## Additional scripts that might be useful

* __Pacbio_get_fasta.pl__   
  _Description_: Extract fasta files (one for ccs.fa, the other for subread.fa) from SMRT runs with two movies (old, pre-2013)
  _Usage_: perl Pacbio_get_fasta.pl <director_of_SMRT_analysis_result> <serial/lane number> <output_filename>
  _Help_: perl Pacbio_get_fasta.pl

* __Pacbio2_get_fasta.pl__   
  _Description_: Extract fasta files (one for ccs.fa, the other for subread.fa) from SMRT runs with two movies (new, post-2013)
  _Usage_: perl Pacbio2_get_fasta.pl <director_of_SMRT_analysis_result> <serial/lane number> <output_filename>
  _Help_: perl Pacbio2_get_fasta.pl
  
* __Pacbio_trim_cDNA_primer2.pl__   
  _Description_: Remove cDNA primers (SMARTer protocol) from the subreads 
  _Usage_: perl Pacbio_trim_cDNA_primer2.pl <subread.fa> <AAGCAGTGGTATCAACGCAGAGTAC> <AAGCAGTGGTATCAACGCAGAGTACATGGG> 0.25 300
  _Help_: perl Pacbio_trim_cDNA_primer2.pl

* __Pacbio_seperate_CCS_2.1.pl__   
  _Description_: Seperate wells with CCS reads with those without CCS reads 
  _Usage_: perl Pacbio_seperate_CCS_2.1.pl <subread.fa.t15> <reads_of_insert.fasta> <SubRead.fa.t15> <CCS.fa> <mySample_1>
  _Help_: perl Pacbio_seperate_CCS_2.pl

* __gtf2refFlat_gmap.pl__   
  _Description_: Convert the gtf format into refFlat format 
  _Usage_: perl gtf2refFlat_gmap.pl <input.gtf> <input.length> <output_refFlat>
  _Help_: perl gtf2refFlat_gmap.pl

* __cluster_refFlat.pl__   
  _Description_: Collapse the PacBio rads and report representative transcript isoforms per gene loci via "clustering" (full-length) pacbio reads, opionally with known exon-junctions from refseq and ensembl annotation, or supporting CAGE data, or requiring at least 2 reads from illumina to support the exon-junctions. The first two arguments are mandatory, the rest can be set to "no" if not available. Note the arguments have to be provided in order as below.
  _Usage_: perl cluster_refFlat.pl <input_refFlat> <output_refFlat.basename> <seqlenth> <strand_info> <coverage> <annotated_junction_bed6_bed12> <Illumina_junction_bed6_bed12> <CAGE_bed6> <min_junction_support> <min_intron_len> <max_intron_len> <min_isof_len> <force_illumina_support>
  _Help_: perl cluster_refFlat.pl


## Tutorial  
1. Prepare Fasta file for PacBio reads
It should be a FASTA file, and we refer to "PacBio.fa" in this tutorial.

2. Prepare Fasta files for Illumina reads and create a file list
As the throughput of NGS grows, using ALL reads in one lane for correction can lead to memory crisis. We can split one big Illumina read file into several smaller ones, perform the correction with the smaller file, and later merge the voting information for a final version.
Let say, we have split the big Illumina file into 26 smaller ones: <NGSa>, <NGSb>, ..., <NGSz>.
We can create a file list with the following command:
```ls -d $PWD/NGS? > NGS_file_list.txt```

3. Prepare parameters
The <SPEC.txt> file is a tab-delimited file, with the first column for keywords recognized by ipec tools, and the second column for user defined files/parameters.
We recommend to generate this file with the following perl script:
```
perl -e 'print join("\t","PacBio_fasta","/Users/arthur/Desktop/IPEC/data/PacBio.fa"),"\n";' > SPEC.txt
perl -e 'print join("\t","Illumina_fasta_list","/Users/arthur/Desktop/IPEC/data/NGS_file_list.txt"),"\n";' >> SPEC.txt
perl -e 'print join("\t","Bowtie2_folder","/Users/arthur/Desktop/bowtie2/"),"\n";' >> SPEC.txt
perl -e 'print join("\t","IPEC_folder","/Users/arthur/Desktop/IPEC/scripts/"),"\n";' >> SPEC.txt
perl -e 'print join("\t","Output_folder","output"),"\n";' >> SPEC.txt
perl -e 'print join("\t","Bowtie2_threads","2"),"\n";' >> SPEC.txt
perl -e 'print join("\t","Remove_tmp_files","yes"),"\n";' >> SPEC.txt
```

4. Prepare pipeline bash file
```perl /PATH_TO_IPEC_FOLDER/scripts/ipec_MAKE.pl SPEC.txt run.sh```

5. Run pipeline
```bash run.sh```

6. RESULT
<output/Results/corrected.fa> is the corrected fasta
Comparing the raw and post-corrected sequences using BLAT to rn6.0:
```
browser details PacBio.fa     2204     5  3111  3113    91.0%  chr2   +    18392251  18416979  24729
browser details corrected.fa  2827     5  2840  2842   100.0%  chr2   +    18392251  18416979  24729
```

7. If you think another iteration is needed, simply use this file as "PacBio_fasta", re-make the <SPEC.txt> file, re-make the "run.sh" and run ipec again.


## Contact
```arthur (dot) yxt (at) gmail (dot) com  ```   for bug reporting or requiring additional functionality


## Citation
1.  Xi Wang, Xintian You, Jingyi Hou, Julian D. Langer, Fiona Rupprecht, Irena Vlatkovic, Claudia Quedenau, Georgi Tushev, Irina Epstein, Bernhard Schaefke, Wei Sun, Liang Fang, Guipeng Li, Yuhui Hu, Erin M Schuman, Wei Chen. (2019) __Full-length transcriptome reconstruction reveals a large diversity of RNA isoforms and open reading frames in the rat hippocampus__, _under review_. 

