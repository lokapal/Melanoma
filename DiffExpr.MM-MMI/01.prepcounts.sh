#!/bin/sh
# Dependencies: featureCounts      https://subread.sourceforge.net/
featureCounts -a /usr/local/genomes/hg38.gtf -o counts.txt -T 30 -t exon -g gene_id \
../MP/rnaseq.MP.rep1.bam ../MP/rnaseq.MP.rep2.bam \
../MM/rnaseq.MM.rep1.bam ../MM/rnaseq.MM.rep2.bam
