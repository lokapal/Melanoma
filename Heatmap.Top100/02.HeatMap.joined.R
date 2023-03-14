#!/usr/bin/Rscript
# script to draw heatmap for two DESeq joined results, selected by top100 up/down regulated genes in MM-MMI case
# (C) Yuri Kravatsky, lokapal@gmail.com, jiri@eimb.ru
# Input:  MP-MM.joined text file, precomputed by diffexp_select_by_list.pl
# Output: 1. MMI-MM-MP.joined.pdf       Complex Heatmap with various clustering approaches.
#
# Dependency tools & libraries:
# 1. R
# Bioconductor packages:
# 2. Complex Heatmap           http://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html
# Common R libraries:
# 3. dendextend, circlize


suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(dendextend))

#GeneID	log2FC	pvalue	padj	Gene
 
data100 <- read.table("MP-MM.joined", sep = "\t", header=T, colClasses = c("character", "numeric", "numeric", "numeric", "character", "numeric", "numeric", "numeric"), row.names=NULL)

dimarr <- nrow(data100)
hei <- dimarr*6/70
cairo_pdf(filename = "MM-MMI.joined.pdf",
width = 3, height = hei, pointsize = 12,
family = "sans", bg = "white", antialias = "default", fallback_resolution = 300, onefile=T)
#png("test2.png", 300, 1200, pointsize=20)

###!!!!!!!heatmap.2(cbind(data100$Exp, data100$Exp), dendrogram="none",labCol="",labRow=data100$Gene,trace="none")!!!!! 
###WORKS!!!!

#mat<-cbind(data100$T,data100$NK)
mat <- cbind(data100$log2FC1, data100$log2FC2)
rownames(mat) <- data100$Gene
colnames(mat) <- c("MM-MMI","MP-MM")
row_dend = as.dendrogram(hclust(dist(mat)))
#row_dend = as.dendrogram(hclust(dist(mat),method = "complete")) ward.D2, mcquitty, median, centroid
row_dend = color_branches(row_dend, k = 4)
dend_thin <- row_dend %>%
     set ("branches_lwd",0.75)

row_dend_complete <- as.dendrogram(hclust(dist(mat),method = "complete")) # ward.D2, mcquitty, median, centroid
row_dend_complete <- color_branches(row_dend_complete, k = 4)
dend_thin_complete <- row_dend_complete %>%
     set ("branches_lwd",0.75)

row_dend_ward <- as.dendrogram(hclust(dist(mat),method = "ward.D2")) # ward.D2, mcquitty, median, centroid
row_dend_ward <- color_branches(row_dend_ward, k = 4)
dend_thin_ward <- row_dend_ward %>%
     set ("branches_lwd",0.75)

row_dend_mcquitty <- as.dendrogram(hclust(dist(mat),method = "mcquitty")) # ward.D2, mcquitty, median, centroid
#row_dend_mcquitty <- color_branches(row_dend_mcquitty, k = 4)
dend_thin_mcquitty <- row_dend_mcquitty %>%
     set ("branches_lwd",0.75)

row_dend_median <- as.dendrogram(hclust(dist(mat),method = "median")) # ward.D2, mcquitty, median, centroid
#row_dend_median <- color_branches(row_dend_median, k = 4)
dend_thin_median <- row_dend_median %>%
     set ("branches_lwd",0.75)

row_dend_centroid <- as.dendrogram(hclust(dist(mat),method = "centroid")) # ward.D2, mcquitty, median, centroid
#row_dend_centroid <- color_branches(row_dend_centroid, k = 4)
dend_thin_centroid <- row_dend_centroid %>%
     set ("branches_lwd",0.75)

#    top_annotation = HeatmapAnnotation(
#        foo = 1:10,
#        annotation_legend_param = list(foo = list(title = "foo_top_anno"))
#    ))

f1 = colorRamp2(c(-5, -2, -1.5, 0, 1.5, 2, 5.5), c("darkblue", "blue", "lightblue", "white", "yellow", "red", "darkred"))

hc <- hclust(dist(mat))
group <- cutree(hc, k = 20)

Heatmap (mat,
row_names_side="left",
heatmap_width = unit(4, "cm"),
cluster_rows = FALSE,
show_column_dend = FALSE,
row_dend_side = "right",
cluster_columns = FALSE,
col=f1,
show_column_names = TRUE,
row_names_gp = gpar(fontsize = 6),
column_names_gp = gpar(fontsize = 6),
heatmap_legend_param = list(
    title = "Log2FC", at = c(-5, -2, -1.5, 0, 1.5, 2, 5.5),
    labels = c("-5.0", "-2.0", "-1.5", "0", "1.5", "2.0", "5.5")
    ),
) 

mat.sorted<-mat[order(mat[,2], decreasing = TRUE), ]

Heatmap (mat.sorted,
row_names_side="left",
heatmap_width = unit(4, "cm"),
cluster_rows = FALSE,
show_column_dend = FALSE,
cluster_columns = FALSE,
col=f1,
show_column_names = TRUE,
row_names_gp = gpar(fontsize = 6),
column_names_gp = gpar(fontsize = 6),
heatmap_legend_param = list(
    title = "Log2FC", at = c(-5, -2, -1.5, 0, 1.5, 2, 5.5),
    labels = c("-5.0", "-2.0", "-1.5", "0", "1.5", "2.0", "5.5")
    ),
) 

Heatmap (mat,
row_names_side="left",
heatmap_width = unit(5, "cm"),
#cluster_rows = TRUE,
#cluster_rows = row_dend,
cluster_rows = dend_thin,
#cluster_rows = cluster_within_group(t(mat), group), 
row_split = 4,
#border = TRUE,
border_gp = gpar(col = "black", lwd = 0.5), ## new version reqd!
show_column_dend = FALSE,
row_dend_side = "right",
cluster_columns = FALSE,
col=f1,
row_names_gp = gpar(fontsize = 6),
column_names_gp = gpar(fontsize = 6),
column_title = "Default clustering",
heatmap_legend_param = list(
    title = "Log2FC", at = c(-5, -2, -1.5, 0, 1.5, 2, 5.5),
    labels = c("-5.0", "-2.0", "-1.5", "0", "1.5", "2.0", "5.5")
    ),
) 


Heatmap (mat,
row_names_side="left",
heatmap_width = unit(5, "cm"),
cluster_rows = dend_thin_complete,
row_split = 4,
border_gp = gpar(col = "black", lwd = 0.5), ## new version reqd!
show_column_dend = FALSE,
row_dend_side = "right",
cluster_columns = FALSE,
col=f1,
row_names_gp = gpar(fontsize = 6),
column_names_gp = gpar(fontsize = 6),
column_title = "Complete clustering",
heatmap_legend_param = list(
    title = "Log2FC", at = c(-5, -2, -1.5, 0, 1.5, 2, 5.5),
    labels = c("-5.0", "-2.0", "-1.5", "0", "1.5", "2.0", "5.5")
    ),
) 

Heatmap (mat,
row_names_side="left",
heatmap_width = unit(5, "cm"),
cluster_rows = dend_thin_ward,
row_split = 4,
border_gp = gpar(col = "black", lwd = 0.5), ## new version reqd!
show_column_dend = FALSE,
row_dend_side = "right",
cluster_columns = FALSE,
col=f1,
row_names_gp = gpar(fontsize = 6),
column_names_gp = gpar(fontsize = 6),
column_title = "Ward2 clustering",
heatmap_legend_param = list(
    title = "Log2FC", at = c(-5, -2, -1.5, 0, 1.5, 2, 5.5),
    labels = c("-5.0", "-2.0", "-1.5", "0", "1.5", "2.0", "5.5")
    ),
) 

Heatmap (mat,
row_names_side="left",
heatmap_width = unit(5, "cm"),
cluster_rows = dend_thin_mcquitty,
#row_split = 4,
#border_gp = gpar(col = "black", lwd = 0.5), ## new version reqd!
show_column_dend = FALSE,
row_dend_side = "right",
cluster_columns = FALSE,
col=f1,
row_names_gp = gpar(fontsize = 6),
column_names_gp = gpar(fontsize = 6),
column_title = "McQuitty clustering",
heatmap_legend_param = list(
    title = "Log2FC", at = c(-5, -2, -1.5, 0, 1.5, 2, 5.5),
    labels = c("-5.0", "-2.0", "-1.5", "0", "1.5", "2.0", "5.5")
    ),
) 

Heatmap (mat,
row_names_side="left",
heatmap_width = unit(5, "cm"),
cluster_rows = dend_thin_median,
#row_split = 4,
#border_gp = gpar(col = "black", lwd = 0.5), ## new version reqd!
show_column_dend = FALSE,
row_dend_side = "right",
cluster_columns = FALSE,
col=f1,
row_names_gp = gpar(fontsize = 6),
column_names_gp = gpar(fontsize = 6),
column_title = "Median clustering",
heatmap_legend_param = list(
    title = "Log2FC", at = c(-5, -2, -1.5, 0, 1.5, 2, 5.5),
    labels = c("-5.0", "-2.0", "-1.5", "0", "1.5", "2.0", "5.5")
    ),
) 

Heatmap (mat,
row_names_side="left",
heatmap_width = unit(5, "cm"),
cluster_rows = dend_thin_centroid,
#row_split = 4,
#border_gp = gpar(col = "black", lwd = 0.5), ## new version reqd!
show_column_dend = FALSE,
row_dend_side = "right",
cluster_columns = FALSE,
col=f1,
row_names_gp = gpar(fontsize = 6),
column_names_gp = gpar(fontsize = 6),
column_title = "Centroid clustering",
heatmap_legend_param = list(
    title = "Log2FC", at = c(-5, -2, -1.5, 0, 1.5, 2, 5.5),
    labels = c("-5.0", "-2.0", "-1.5", "0", "1.5", "2.0", "5.5")
    ),
) 
 
dev.off();
