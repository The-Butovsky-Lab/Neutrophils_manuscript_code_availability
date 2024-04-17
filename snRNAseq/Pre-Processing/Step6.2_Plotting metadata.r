# Donors have been genotyped for APOE, and have tau and amyloid pathology scores from the same tissue blocks used for snRNAseq. CellRanger count martices
# are stored in the input folder, with a seperate folder for each donor. The metadata .csv file is in the input folder. The output folder is created. 

# 5.1. Plot metadata of the dataset

# Author: Thomas Rust (THR) adapted from Mirjam Koster (MKO)
# Date: 26/10/2023 (updated)

# load packages
library(Seurat)
library(ggplot2)
library(scales)
library(dplyr)

# set output directory and load in seurat object
Output.dir <- "Output/"
total <- readRDS(file = paste0(Output.dir, "5.0 Total object full metadata.rds"))

# set subdirectory for metadata plots
Output.dir.meta <- paste0(Output.dir, "5.1 Plotting metadata/")

## nCount and nFeature need to be log transformed

total$nCount_RNA_log <- log(total$nCount_RNA)
total$nFeature_RNA_log <- log(total$nFeature_RNA)
to_plot <- c("nCount_RNA_log", "nCount_RNA", "nFeature_RNA", "nFeature_RNA_log")

for(i in to_plot){
  png(filename = paste0(Output.dir.meta, "5. ",i,".png"),  width = 1080, height = 1080)
print(
    FeaturePlot(object = total, features = i, raster = FALSE, label = TRUE, min.cutoff = "q1", max.cutoff = "q99",
            cols = c("grey90", "#009DFF", "#C929DA","#D40000", "darkred"))
  )
dev.off()
}

png(filename = paste0(Output.dir.meta, "5. nCount_log nFeature.png"),  width = 1920, height = 1080)

  FeaturePlot(object = total, features = "nCount_RNA_log", raster = FALSE, label = TRUE, min.cutoff = "q1", max.cutoff = "q99",
            cols = alpha(c("grey90", "#009DFF", "#C929DA","#D40000", "darkred"),0.5), pt.size = 1/1000000) -
    FeaturePlot(object = total, features = "nFeature_RNA", raster = FALSE, label = TRUE, min.cutoff = "q1", max.cutoff = "q99",
                cols = alpha(c("grey90", "#009DFF", "#C929DA","#D40000", "darkred"),0.5), pt.size = 1/1000000)
dev.off()

#Categorical variables
total@meta.data$apoe <- as.character(total@meta.data$apoe)
total@meta.data$RIN <- as.numeric(total@meta.data$RIN)
to_plot_cat <- c("group", "diagnosis", "apoe", "sex", "snRNAseq_batch", "pmd_min")

for(i in to_plot_cat){
  png(filename = paste0(Output.dir.meta, "5. Metadata plot: ",i,".png"),  width = 1080, height = 1080)
  j = length(unique(total@meta.data[[i]]))
  print(
    DimPlot(object = total, group.by = i, raster = FALSE, 
            cols = alpha(hue_pal()(j), 0.5))
  )
dev.off()
}

# Continuous variables
to_plot_con <- c("age", "RIN", "percent.ribo", "percent.mito", "pmd_min")

for(i in to_plot_con){
  png(filename = paste0(Output.dir.meta, "5. Metadata plot continuous: ",i,".png"),  width = 1080, height = 1080)
  print(
    FeaturePlot(object = total, features = i, raster = FALSE, min.cutoff = "q1", max.cutoff = "q99",
            cols = alpha(c("grey90", "#009DFF", "#C929DA","#D40000", "darkred"),0.5))
  )
  dev.off()
}

# Proceed to THR 6.0 Subsetting all cell types