# Donors have been genotyped for APOE, and have tau and amyloid pathology scores from the same tissue blocks used for snRNAseq. CellRanger count martices
# are stored in the input folder, with a seperate folder for each donor. The metadata .csv file is in the input folder. The output folder is created. 

# 5.0. Update seurat object metadata

# Author: Thomas Rust (THR) 
# Date: 23/02/2024 (updated)

# load packages
library(ggplot2)
library(Seurat)
library(Matrix)
library(reshape2)
library(dplyr)

#set output directory and load in total object
Output.dir <- "Output/"
total <- readRDS(file = paste0(Output.dir, "4.2 all_samples_cell_category.rds"))
total
total@meta.data$Cell_category

# Add updated metadata with pmd scores 

metadata0 <- read.csv2("Input/metadata_apoe33_34_pmd.csv")
metadata0$sample <- as.character(metadata0$sample)
metadata0$pmd_min <- gsub(",", ".", metadata0$pmd_min)
metadata0$pmd_min <- as.numeric(metadata0$pmd_min)
metadata0

# Add the metadata to each individual sample in the list                 
metadata1 <- left_join(total[["sample"]], metadata0)     
row.names(metadata1) <- row.names(total[[]])
total <- AddMetaData(total, metadata = metadata1)
total@meta.data

# save final object for analysis
saveRDS(object = total, file = paste0(Output.dir, "5.0 Total object full metadata.rds"))