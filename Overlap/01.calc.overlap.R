#!/usr/bin/Rscript
# script to calculate the number of interaction of many lists by SuperExactTest:
# Wang, M.; Zhao, Y.; Zhang, B. Efficient test and visualization of multi-set intersections. Sci Rep 2015, 5, 16923.
# (C) Yuri Kravatsky, lokapal@gmail.com, jiri@eimb.ru
# Input:  list.txt with the following columns:
#    1. DSB HEK293 top 4921 genes from DOI: 10.3390/ijms23137201
#    2. 4C  HEK293 top 4921 genes from DOI: 10.3390/cells8111393
#    3. Mel Z top 101 upregulated   genes by Log2FC > 1.5 
#    4. Mel Z top 122 downregulated genes by Log2FC < -1.5
#    5. hg38 EMBL v106 Gene IDs complete list
#    6. hg38 EMBL v106 Gene Names complete list
# Output: 1. summary.table.csv    Complete results, including lists of overlapping genes and p-values
#
# Dependency tools & libraries:
# 1. R
# CRAN packages:
# 2. SuperExactTest 

# read lists.txt file, convert empty cells to NAs
lists <- read.table("lists.txt", sep = "\t", header=T, na.strings=c("", "NA"), row.names=NULL)

# Load all lists from variables, omitting NAs
DSBs    <- as.character(na.omit(lists$DSB))
FCs     <- as.character(na.omit(lists$FourC))
upreg   <- as.character(na.omit(lists$UPREG))
downreg <- as.character(na.omit(lists$DOWNREG))
# Get "universe" size by calculating genome gene list size
hgg     <- length(as.character(na.omit(lists$GeneName)))

mylist <- list(DSBs,FCs, upreg, downreg)
names(mylist) <- c("DSB", "FourC", "UpReg", "DownReg")

suppressPackageStartupMessages(library('SuperExactTest'))
res=supertest(mylist, n=hgg)
summary(res)
plot(res, Layout="landscape", degree=2:4, sort.by="p-value", margin=c(0.5,5,1,2), log.scale=TRUE, min.intersection.size=1)
write.csv(summary(res)$Table, file="summary.table.csv", row.names=FALSE)
