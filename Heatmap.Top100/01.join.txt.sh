#!/bin/sh
# script to create source data for further heatmap drawing
# the purpose - to select DESeq2 results from MP-MM.tsv
# according to the gene IDs that are present in MM-MMI.top100.txt
# output: file MP-MM.joined text file
perl diffexp_select_by_list.pl MP-MM.tsv MM-MMI.top100.txt
