#!/bin/sh
# script to get and process H.sapiens Mel Z grown on Matrigel (MM) hg38 gene expression in TPM
# (C) Yuri Kravatsky, lokapal@gmail.com, jiri@eimb.ru
# Input:  H.sapiens Mel Z cells grown on Matrigel (MM) RNA-Seq data: GEO GSE221872: SRR18867216, SRR18867217
# Output: 1. MM.RNAseq1.genes.results    Mel Z MM hg38 gene expression replica 1
#            MM.RNAseq2.genes.results    Mel Z MM hg38 gene expression replica 1
#         2. MM.hg38.TPM                 joined replicas expression all and average values
#            MM.hg38.mean.TPM            Mel Z MM hg38 gene expression average values per gene
#
# Dependency tools:
# 1. NCBI SRA Toolkit     https://github.com/ncbi/sra-tools
# 2. Trimmomatic 0.39     http://www.usadellab.org/cms/?page=trimmomatic
# 3. RSEM 1.3.1           https://github.com/deweylab/RSEM
# 4. STAR 2.7.5c          https://github.com/alexdobin/STAR
# 5. R with libraries     tibble, dplyr
# Get data from NCBI SRA Archive
fastq-dump --split-files --gzip SRR18867216
fastq-dump --split-files --gzip SRR18867217
mv SRR18867216_1.fastq.gz MM.rep1.fastq.gz
mv SRR18867217_1.fastq.gz MM.rep2.fastq.gz
# Quality trimming
trimmomatic SE -threads 30 MM.rep1.fastq.gz rep1.filtered.fastq.gz LEADING:18 TRAILING:18 SLIDINGWINDOW:4:22 MINLEN:20
trimmomatic SE -threads 30 MM.rep2.fastq.gz rep2.filtered.fastq.gz LEADING:18 TRAILING:18 SLIDINGWINDOW:4:22 MINLEN:20
# Expression calculation per replicate
rsem-calculate-expression -p 32 --fragment-length-mean 255 --star -p 32 --output-genome-bam --calc-ci \
--ci-memory 30720 --star-gzipped-read-file rep1.filtered.fastq.gz /usr/local/genomes/hg38.rsem/hg38 MM.RNAseq1
rsem-calculate-expression -p 32 --fragment-length-mean 255 --star -p 32 --output-genome-bam --calc-ci \
--ci-memory 30720 --star-gzipped-read-file rep2.filtered.fastq.gz /usr/local/genomes/hg38.rsem/hg38 MM.RNAseq2
# Calculate average expression values between replicas
./lib/mean_expr.R
rm -f *.bam
