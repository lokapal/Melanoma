#!/usr/bin/Rscript

suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))
namehead = "MM.RNAseq"
nametail = ".genes.results"

for (num in 1:2) {
   fname <- paste0(namehead,num,nametail)
#  cat(fname,"\n")
   temptbl <- read.table(fname, skip=1, sep = "\t", strip.white = TRUE)
   newname <- paste0("Expr",num)
   FPKM    <- paste0("FPKM",num)
   colnames(temptbl)[6] <- newname
   colnames(temptbl)[7] <- FPKM
   if (num == 1) {
         alldata <- as.data.frame (temptbl$newname,row.names=temptbl$V1)
         alldata <- as.data.frame (temptbl[[newname]],row.names=as.vector(temptbl$V1))
         colnames(alldata)[num] <- newname
#         head(alldata)

         FPdata <- as.data.frame (temptbl$FPKM,row.names=temptbl$V1)
         FPdata <- as.data.frame (temptbl[[FPKM]],row.names=as.vector(temptbl$V1))
         colnames(FPdata)[num] <- FPKM
         next
                 }
   alldata <- add_column(alldata, temptbl[[newname]])
   FPdata  <- add_column(FPdata,  temptbl[[FPKM]])
   colnames(alldata)[num] <- newname
   colnames(FPdata)[num]  <- FPKM
                 }

  medvect  <- apply(alldata, 1, median)
#  meanvect <- apply(alldata, 1, mean, trim=0.2) # trimmed mean
  meanvect <- apply(alldata, 1, mean)          # "true" mean


  meanpkm  <- apply(FPdata, 1, mean)          # "true" mean
  medpkm   <- apply(FPdata, 1, median)          # "true" mean

#  collist <- c(2:9)
#  subset <- alldata %>% select(collist)
#write.table(as.data.frame(subset), file="subset.txt", row.names=F, col.names=T, sep="\t", quote=F)

  alldata <- add_column(alldata,medvect,meanvect)
  alldata <- rownames_to_column(alldata, var = "GeneID") %>% as_tibble()
  finalmed  <- alldata %>% select(GeneID,medvect)
  finalmean <- alldata %>% select(GeneID,meanvect)

  FPdata <-  add_column(FPdata,medpkm,meanpkm)
  FPdata <-  rownames_to_column(FPdata, var = "GeneID") %>% as_tibble()
  finalPKmed  <- FPdata %>% select(GeneID,medpkm)
  finalPKmean <- FPdata %>% select(GeneID,meanpkm)

colnames(finalPKmean)[2] <- "FPKM, mean"
colnames(finalmean)[2]   <- "TPM, mean"

write.table(as.data.frame(alldata), file="MM.hg38.TPM", row.names=F, col.names=T, sep="\t", quote=F)
write.table(as.data.frame(FPdata), file="MM.hg38.FPKM", row.names=F, col.names=T, sep="\t", quote=F)
write.table(as.data.frame(finalmean), file="MM.hg38.mean.TPM", row.names=F, col.names=T, sep="\t", quote=F)
write.table(as.data.frame(finalPKmean), file="MM.hg38.mean.FPKM", row.names=F, col.names=T, sep="\t", quote=F)
