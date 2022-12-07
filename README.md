# Yao_RNA_seq
## 1.1 Data preparation

**Option 1**: You can download from the link they provided.
For Ruddle users, the following line will create symlinks to your samples:
```
/home/bioinfo/software/knightlab/bin_Mar2018/ycgaFastq  fcb.ycga.yale.edu:3010/HqlP2H0hKOLOqu0y4PJZawtU6UwEk/sample_dir_000009572
```
**Option 2**: Download your raw data to your laptop and then transfer them to the web sever via Globus.

## 1.2 Data rename and unzip
```
nohup gunzip -c Rse1-1_164_029_S59_L002_R1_001.fastq.gz > G1_Rep1_R1 &
nohup gunzip -c Rse1-1_164_029_S59_L002_R2_001.fastq.gz > G1_Rep1_R2 &
nohup gunzip -c Rse1-2_152_041_S60_L002_R1_001.fastq.gz > G1_Rep2_R1 &
nohup gunzip -c Rse1-2_152_041_S60_L002_R2_001.fastq.gz > G1_Rep2_R2 &
nohup gunzip -c Rse1-3_140_053_S61_L002_R1_001.fastq.gz > G1_Rep3_R1 &
nohup gunzip -c Rse1-3_140_053_S61_L002_R2_001.fastq.gz > G1_Rep3_R2 &

nohup gunzip -c Rse1-4_128_065_S62_L002_R1_001.fastq.gz > G2_Rep1_R1 &
nohup gunzip -c Rse1-4_128_065_S62_L002_R2_001.fastq.gz > G2_Rep1_R2 &
nohup gunzip -c Rse1-5_116_077_S63_L002_R1_001.fastq.gz > G2_Rep2_R1 &
nohup gunzip -c Rse1-5_116_077_S63_L002_R2_001.fastq.gz > G2_Rep2_R2 &
nohup gunzip -c Rse1-6_104_089_S64_L002_R1_001.fastq.gz > G2_Rep3_R1 &
nohup gunzip -c Rse1-6_104_089_S64_L002_R2_001.fastq.gz > G2_Rep3_R2 &
```

## 2. QC and Data trimming
Generally, we need to check the quailty of our raw data and trim the low-quaility data.
### 2.1 QC
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
For the interpretation of QC reports, you can refer to the following materials.

Chinese version: https://blog.csdn.net/weixin_43569478/article/details/108079243

English version: https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/qc_fastqc_assessment.html

### 2.2 Trimming
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
### 3.1 Preparing the human reference genome
https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.40
```
nohup gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz &
nohup gunzip GCF_000001405.40_GRCh38.p14_genomic.gtf.gz &
```

### 3.2 Mapping to the reference genome
```

vi run_paired_hisat2_mapping_self_indexed_RefGenome.sh

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

## 4. Reads counting
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

