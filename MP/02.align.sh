#!/bin/sh
# script to align RNASeq data with the same parameters as in rsem 1.3.1
# (C) Yuri Kravatsky, lokapal@gmail.com, jiri@eimb.ru
# Input:  pre-filtered RNASeq rep1 and RNASeq rep2 datasets
# Output: 1. rnaseq.MP.rep1.bam          Mel Z MP RNASeq replica 1 alignment to hg38
#            rnaseq.MP.rep2.bam          Mel Z MP RNASeq replica 2 alignment to hg38
#         2. MP_Pearson.pdf              Mel Z MP RNASeq replica consistency check: Pearson correlation
#            MP_Spearman.pdf             Mel Z MP RNASeq replica consistency check: Spearman correlation
#
# Dependency tools:
# 1. STAR 2.7.5c          https://github.com/alexdobin/STAR
# 2. samtools 1.15        http://www.htslib.org
# 3. deepTools2 3.5.1     https://deeptools.readthedocs.io
# align to hg38 genome with the options that are the same as in rsem 1.3.1
STAR --genomeDir /usr/local/genomes/hg38.rsem --outSAMunmapped Within --outFilterType BySJout \
--outSAMattributes NH HI AS NM MD --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
--alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --runThreadN 32 \
--genomeLoad NoSharedMemory --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM \
--outFileNamePrefix rep1. --readFilesCommand zcat --readFilesIn rep1.filtered.fastq.gz
# 
STAR --genomeDir /usr/local/genomes/hg38.rsem --outSAMunmapped Within --outFilterType BySJout \
--outSAMattributes NH HI AS NM MD --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
--alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --runThreadN 32 \
--genomeLoad NoSharedMemory --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM \
--outFileNamePrefix rep2. --readFilesCommand zcat --readFilesIn rep2.filtered.fastq.gz
mv rep1.Aligned.sortedByCoord.out.bam rnaseq.MP.rep1.bam
samtools index -@ 16 rnaseq.MP.rep1.bam
mv rep2.Aligned.sortedByCoord.out.bam rnaseq.MP.rep2.bam
samtools index -@ 16 rnaseq.MP.rep2.bam
# Quality Control to check consistency between replicas (optional)
multiBamSummary bins -p 30 -b rnaseq.MP.rep1.bam rnaseq.MP.rep2.bam -o res12.npz
plotCorrelation -in res12.npz --skipZeros -c spearman --removeOutliers -p scatterplot -o MP_Spearman.pdf --labels rep1 rep2 --log1p
plotCorrelation -in res12.npz --skipZeros -c pearson  --removeOutliers -p scatterplot -o MP_Pearson.pdf  --labels rep1 rep2 --log1p
