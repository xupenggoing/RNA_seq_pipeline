# RNA_seq_pipeline
This is the pipeline for our RNA-seq data analysis. There are two groups, G1 and G2; for each group, there are three replicates, Rep1, Rep2, and Rep3; each replicate contains R1 and R2 reads. The command lines are preliminary. Of course, both the command lines and file-naming could be better.
## 1. Data preparation
### 1.1 Load datasets to web server

- **Option 1**: You can download from the link the sequence service provider gave to youp.
- **Option 2**: Download your raw data to your laptop and then transfer them to the web sever via Globus.

### 1.2 Rename and unzip datasets
Generally, the names for your raw data are too long and complicated to read and for the following analysis. We'd better change their names. We can finish unzipping and name-changing in a single step. Note here `-c` means keeping the raw data; otherwise, the raw data will dispear and only the name-changed files will be kept.
```
nohup gunzip -c long_file_name_R1_001.fastq.gz > G1_Rep1_R1 &
nohup gunzip -c long_file_name_R2_001.fastq.gz > G1_Rep1_R2 &
nohup gunzip -c long_file_name_R1_001.fastq.gz > G1_Rep2_R1 &
nohup gunzip -c long_file_name_R2_001.fastq.gz > G1_Rep2_R2 &
nohup gunzip -c long_file_name_R1_001.fastq.gz > G1_Rep3_R1 &
nohup gunzip -c long_file_name_R2_001.fastq.gz > G1_Rep3_R2 &

nohup gunzip -c long_file_name_R1_001.fastq.gz > G2_Rep1_R1 &
nohup gunzip -c long_file_name_R2_001.fastq.gz > G2_Rep1_R2 &
nohup gunzip -c long_file_name_R1_001.fastq.gz > G2_Rep2_R1 &
nohup gunzip -c long_file_name_R2_001.fastq.gz > G2_Rep2_R2 &
nohup gunzip -c long_file_name_R1_001.fastq.gz > G2_Rep3_R1 &
nohup gunzip -c long_file_name_R2_001.fastq.gz > G2_Rep3_R2 &
```
## 2. QC and Data trimming
Generally, we need to check the quailty of our raw data and trim the low-quaility data. 
### 2.1 QC
FastQC, a java written software, is developed by the Bioinformatics Group at the Babraham Institute, UK. It's thought to be the most used software for next generation sequences QC.

For yale web server users, firstly you need to load the software you need for analysis via `module load FastQC/0.11.9-Java-11`. However, you need to know which software we have in the yale web server via `module avail` to see all the available softwares.
```
module load FastQC/0.11.9-Java-11
nohup fastqc G1_Rep1_R1.fastq &
nohup fastqc G1_Rep1_R2.fastq &
nohup fastqc G1_Rep2_R1.fastq &
nohup fastqc G1_Rep2_R2.fastq &
nohup fastqc G1_Rep3_R1.fastq &
nohup fastqc G1_Rep3_R2.fastq &

nohup fastqc G2_Rep1_R1.fastq &
nohup fastqc G2_Rep1_R2.fastq &
nohup fastqc G2_Rep2_R1.fastq &
nohup fastqc G2_Rep2_R2.fastq &
nohup fastqc G2_Rep3_R1.fastq &
nohup fastqc G2_Rep3_R2.fastq &
```
After the QC, we will get two files, one ended up with `.html` and the other with `.zip`. For `.html` files containing the QC report, you can transfer them to your laptop and open them with your web browser. For the interpretation of QC reports, you can refer to: https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/qc_fastqc_assessment.html

### 2.2 Trimming
When trimming your data, you need to figure out whether your sequencing reads are pair-end, it's easy to confrim. We used Illumina sequencing machine, thus I put `--illumina` here.
```
module load Trim_Galore/0.6.6-GCCcore-10.2.0-Python-3.8.6
nohup trim_galore --paired --illumina G1_Rep1_R1.fastq G1_Rep1_R2.fastq &
nohup trim_galore --paired --illumina G1_Rep2_R1.fastq G1_Rep2_R2.fastq &
nohup trim_galore --paired --illumina G1_Rep3_R1.fastq G1_Rep3_R2.fastq &

nohup trim_galore --paired --illumina G2_Rep1_R1.fastq G2_Rep1_R2.fastq &
nohup trim_galore --paired --illumina G2_Rep2_R1.fastq G2_Rep2_R2.fastq &
nohup trim_galore --paired --illumina G2_Rep3_R1.fastq G2_Rep3_R2.fastq &
```
## 3. Mapping
Mapping is the most important step for RNA-seq analysis.
### 3.1 Preparing the human reference genome
You can download the human reference genome from the following link: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.40. Note that you need to download both reference genome data (.fna) and annotation data (.gft). You can download them to you laptop and then tranfer then to the web server. Also, you can download them to the web server directly via `wget link_address`.

```
wget https://www.ncbi.nlm.nih.gov/projects/r_gencoll/ftp_service/nph-gc-ftp-service.cgi/?HistoryId=MCID_6390e4bdef461250ef424de0&QueryKey=1&ReleaseType=RefSeq&FileType=GENOME_FASTA&Flat=true

wget https://www.ncbi.nlm.nih.gov/projects/r_gencoll/ftp_service/nph-gc-ftp-service.cgi/?HistoryId=MCID_6390e4bdef461250ef424de0&QueryKey=1&ReleaseType=RefSeq&FileType=GENOME_GTF&Flat=true

nohup gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz &
nohup gunzip GCF_000001405.40_GRCh38.p14_genomic.gtf.gz &
```
Then you need to build the reference genome because they can't be used directly.
```
nohup hisat2-build GCF_000001405.40_GRCh38.p14_genomic.fna hisat2_built_Genome 1>hisat2-build.log 2>&1 &
```
After the building, you will get 8 files.
### 3.2 Mapping reads to human reference genome
```
vi run_paired_hisat2_mapping_self_indexed_RefGenome.sh

##Creat a hisat2 mapping script - start line##

#!/bin/bash
#SBATCH --job-name=run_paired_hisat2_mapping_self_indexed_RefGenome
#SBATCH --out="slurm-%j.out"
#SBATCH --time=1-
#SBATCH --nodes=1
#SBATCH -c 20
#SBATCH --mem=64G
#SBATCH --mail-type=ALL

hisat2 --new-summary -p 10 -x ../Human_Ref_Genome/hisat2_built_Genome -1 G1_Rep1_R1_val_1.fq -2 G1_Rep1_R2_val_2.fq -S G1_Rep1_selfIndex.sam --rna-strandness RF --summary-file G1_Rep1_selfIndex_alignment.txt
hisat2 --new-summary -p 10 -x ../Human_Ref_Genome/hisat2_built_Genome -1 G1_Rep2_R1_val_1.fq -2 G1_Rep2_R2_val_2.fq -S G1_Rep2_selfIndex.sam --rna-strandness RF --summary-file G1_Rep2_selfIndex_alignment.txt
hisat2 --new-summary -p 10 -x ../Human_Ref_Genome/hisat2_built_Genome -1 G1_Rep3_R1_val_1.fq -2 G1_Rep3_R2_val_2.fq -S G1_Rep3_selfIndex.sam --rna-strandness RF --summary-file G1_Rep3_selfIndex_alignment.txt

hisat2 --new-summary -p 10 -x ../Human_Ref_Genome/hisat2_built_Genome -1 G2_Rep1_R1_val_1.fq -2 G2_Rep1_R2_val_2.fq -S G2_Rep1_selfIndex.sam --rna-strandness RF --summary-file G2_Rep1_selfIndex_alignment.txt
hisat2 --new-summary -p 10 -x ../Human_Ref_Genome/hisat2_built_Genome -1 G2_Rep2_R1_val_1.fq -2 G2_Rep2_R2_val_2.fq -S G2_Rep2_selfIndex.sam --rna-strandness RF --summary-file G2_Rep2_selfIndex_alignment.txt
hisat2 --new-summary -p 10 -x ../Human_Ref_Genome/hisat2_built_Genome -1 G2_Rep3_R1_val_1.fq -2 G2_Rep3_R2_val_2.fq -S G2_Rep3_selfIndex.sam --rna-strandness RF --summary-file G2_Rep3_selfIndex_alignment.txt

##Creat a hisat2 mapping script - end line##

module load HISAT2/2.2.1-gompi-2020b
sbatch run_paired_hisat2_mapping_self_indexed_RefGenome.sh
```
### 3.3 SAM files to BAM files
```
vi run_selfIndex_sam_to_bam.sh

#!/bin/bash
#SBATCH --job-name=run_selfIndex_sam_to_bam
#SBATCH --out="slurm-%j.out"
#SBATCH --time=1-
#SBATCH --nodes=1
#SBATCH -c 20
#SBATCH --mem=64G
#SBATCH --mail-type=ALL

samtools sort G1_Rep1_selfIndex.sam -o G1_Rep1_selfIndex.bam
samtools sort G1_Rep2_selfIndex.sam -o G1_Rep2_selfIndex.bam
samtools sort G1_Rep3_selfIndex.sam -o G1_Rep3_selfIndex.bam

samtools sort G2_Rep1_selfIndex.sam -o G2_Rep1_selfIndex.bam
samtools sort G2_Rep2_selfIndex.sam -o G2_Rep2_selfIndex.bam
samtools sort G2_Rep3_selfIndex.sam -o G2_Rep3_selfIndex.bam

module load SAMtools/1.16-GCCcore-10.2.0
sbatch run_selfIndex_sam_to_bam.sh
```
### 3.4 Indexing BAM files
```
vi run_indexing_selfIndexedRefGenome_bam.sh

#!/bin/bash
#SBATCH --job-name=run_indexing_selfIndexedRefGenome_bam
#SBATCH --out="slurm-%j.out"
#SBATCH --time=1-
#SBATCH --nodes=1
#SBATCH -c 20
#SBATCH --mem=64G
#SBATCH --mail-type=ALL

samtools index G1_Rep1_selfIndex.bam
samtools index G1_Rep2_selfIndex.bam
samtools index G1_Rep3_selfIndex.bam

samtools index G2_Rep1_selfIndex.bam
samtools index G2_Rep2_selfIndex.bam
samtools index G2_Rep3_selfIndex.bam

module load SAMtools/1.16-GCCcore-10.2.0
sbatch run_indexing_selfIndexedRefGenome_bam.sh 
```
## 4. Reads-counting by featureCounts
```
vi featureCounts_via_Subread_selfIndexRefGenome.sh

#!/bin/bash
#SBATCH --job-name=featureCounts_via_Subread_selfIndexRefGenome
#SBATCH --out="slurm-%j.out"
#SBATCH --time=1-
#SBATCH --nodes=1
#SBATCH -c 20
#SBATCH --mem=64G
#SBATCH --mail-type=ALL

featureCounts -p -T 24 -t exon -g gene_id -a ../Human_Ref_Genome/GCF_000001405.40_GRCh38.p14_genomic.gtf -o Yao_RNA_counts.txt G1_Rep1_selfIndex.bam G1_Rep2_selfIndex.bam G1_Rep3_selfIndex.bam G2_Rep1_selfIndex.bam G2_Rep2_selfIndex.bam G2_Rep3_selfIndex.bam

sbatch featureCounts_via_Subread_selfIndexRefGenome.sh
```

After getting the expression matrix, we can refer to other tutorials, such as http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html.
