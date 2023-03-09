#!/bin/sh
# Dependencies: featureCounts      https://subread.sourceforge.net/
featureCounts -a /usr/local/genomes/hg38.gtf -o counts.txt -T 30 -t exon -g gene_id \
../MM/rnaseq.MM.rep1.bam ../MM/rnaseq.MM.rep2.bam \
../MMI/rnaseq.MMI.rep1.bam ../MMI/rnaseq.MMI.rep2.bam
