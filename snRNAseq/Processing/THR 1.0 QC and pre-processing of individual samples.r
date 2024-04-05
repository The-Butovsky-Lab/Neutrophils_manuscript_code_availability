# Analysis of snRNA-seq data for collaboration with Oleg Butovsky. This analysis conists of only APOE33 and APOE34 donors. Samples are fresh-frozen brain tissue from AD and non-demented CTRL.
# Donors have been genotyped for APOE, and have tau and amyloid pathology scores from the same tissue blocks used for snRNAseq. CellRanger count martices
# are stored in the input folder, with a seperate folder for each donor. The metadata334 .csv file is in the input folder. The output folder is created. 
# 1. QC and processing of individual samples prior to doublet exlcusion and sample integration
# Author: Thomas Rust (THR) adapted from Mirjam Koster (MKO)
# Date: 14/01/2024 (updated)

library(Matrix)
library(Seurat)
library(reticulate)
library(dplyr)
library(patchwork)
library(ggplot2)
sessionInfo()

Output.dir <- "Output/"
total.list <- readRDS(file = paste0(Output.dir, "0. listed_samples unprocessed with metadata.rds"))

#############################################################################
### 1. QC and pre-processing of individual samples
#############################################################################

## QUALITY CONTROL
pdf(file = paste0(Output.dir, "1.1 Counts-Features and histogram per sample.pdf"))
for (a in 1:length(total.list)) {
  sample_name <- names(total.list[a])
  
  plot(x = total.list[[a]]$nCount_RNA, y = total.list[[a]]$nFeature_RNA)
  title(sample_name)
  
  hist(as.matrix(log10(total.list[[a]]@assays$RNA@counts)), main = sample_name)
}
dev.off()

## PRE-PROCESSING
# Percentage mitochondrial & ribosomal genes (metadata)
for(a in 1:length(total.list)) {
  total.list[[a]][["percent.mito"]] <- PercentageFeatureSet(total.list[[a]], pattern = "^MT-")
  total.list[[a]][["percent.ribo"]] <- PercentageFeatureSet(total.list[[a]], pattern = "^RP[SL]")
}

# Visualize QC metrics as a violin plot
pdf(file = paste0(Output.dir,"1.2 Quality control plots.pdf"), 
    title = "Quality control plots")
for(a in 1:length(total.list)) {
  print(VlnPlot(total.list[[a]], features = c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo"), 
                ncol = 2, pt.size = 0))
}
dev.off()


# Subset cells with <5% mito genes and based on n_Feature and n_Count cut-offs
for(a in 1:length(total.list)) {
  total.list[[a]] <- subset(total.list[[a]], subset = percent.mito < 5 & nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 25000 & nCount_RNA > 200) 
}

##NORMALISATION
# default: normalization.method = "LogNormalize", scale.factor = 10000
pdf(file = paste0(Output.dir, "1.3 Normalised data.pdf"))
for(a in 1:length(total.list)) {
  total.list[[a]] <- NormalizeData(total.list[[a]])
  
  hist(as.matrix(log2(total.list[[a]][["RNA"]]@data)), main = paste0("Normalised data: sample ", names(total.list[a])))
}
dev.off()

##FEATURE SELECTION
# Identification of highly variable features
# default: selection.method = "vst", nfeature = 2000

pdf(file = paste0(Output.dir, "1.4 Variable features plot.pdf"))
for(a in 1:length(total.list)) {
  total.list[[a]] <- FindVariableFeatures(total.list[[a]])
  
  # plot variable features with top 10 labeled
  top10 <- head(VariableFeatures(total.list[[a]]), 10)
  plot1 <- VariableFeaturePlot(total.list[[a]])
  
  print(LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0))
}
dev.off()

##SCALE FEATURES
# Scale data & run PCA separately per sample
dim = 1:30 #Check if 30 dimensions is enough by inspecting the Elbow plots

pdf(file = paste0(Output.dir, "1.5 Scaled data + plots per sample.pdf"))
for(a in 1:length(total.list)) {
  total.list[[a]] <- ScaleData(total.list[[a]])
  total.list[[a]] <- RunPCA(total.list[[a]])
  
  print(DimPlot(object = total.list[[a]], reduction = "pca"))
  print(ElbowPlot(total.list[[a]], ndims = 50))
  
  total.list[[a]] <- FindNeighbors(total.list[[a]], dims = dim)
  total.list[[a]] <- FindClusters(object = total.list[[a]], resolution = 0.1)
  total.list[[a]] <- RunTSNE(total.list[[a]], dims = dim)
  total.list[[a]] <- RunUMAP(total.list[[a]], dims = dim)
  
  print(DimPlot(object = total.list[[a]], reduction = "tsne", raster = FALSE))
  print(DimPlot(object = total.list[[a]], reduction = "umap", raster = FALSE))
}
dev.off()

#Save the list of individually-processed samples and proceed to doublet exclusion with the next script (THR 2.0 Doublet exclusion)
rm(list=setdiff(ls(), c("total.list", "Output.dir")))
saveRDS(total.list, file = paste0(Output.dir, "1. listed_samples processed.rds"))