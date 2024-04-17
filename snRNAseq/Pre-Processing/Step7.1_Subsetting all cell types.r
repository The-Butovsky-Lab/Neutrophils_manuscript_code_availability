# Analysis of snRNA-seq data for collaboration with Oleg Butovsky. Samples are fresh-frozen brain tissue from AD and non-demented CTRL.
# Donors have been genotyped for APOE, and have tau and amyloid pathology scores from the same tissue blocks used for snRNAseq. CellRanger count martices
# are stored in the input folder, with a seperate folder for each donor. The metadata .csv file is in the input folder. The output folder is created.

# 6.0 Subsetting of cell types into seperate rds objects. 
# I will then perform integration, sclaing and clustering for each cell type object. 
# Author: Thomas Rust (THR)
# Date: 26/10/2023 (updated)

# load packages
library(ggplot2)
library(Seurat)
library(Matrix)
library(reshape2)
library(dplyr)

# set output directory and read in seurat object for total dataset
Output.dir <- "Output/"
total <- readRDS(file = paste0(Output.dir, "5.0 Total object full metadata.rds"))

# Remove previous clustering information for each resolution as this will need to be recalculated for each subclustering
total@meta.data$integrated_snn_res.0 <- NULL 
total@meta.data$integrated_snn_res.0.01 <- NULL
total@meta.data$integrated_snn_res.0.02 <- NULL 
total@meta.data$integrated_snn_res.0.03 <- NULL
total@meta.data$integrated_snn_res.0.04 <- NULL
total@meta.data$integrated_snn_res.0.05 <- NULL
total@meta.data$integrated_snn_res.0.06 <- NULL 
total@meta.data$integrated_snn_res.0.07 <- NULL 
total@meta.data$integrated_snn_res.0.08 <- NULL 
total@meta.data$integrated_snn_res.0.09 <- NULL 
total@meta.data$integrated_snn_res.0.1 <- NULL 
total@meta.data$integrated_snn_res.0.11 <- NULL
total@meta.data$integrated_snn_res.0.12 <- NULL 
total@meta.data$integrated_snn_res.0.13 <- NULL 
total@meta.data$integrated_snn_res.0.14 <- NULL 
total@meta.data$integrated_snn_res.0.15 <- NULL 
total@meta.data$integrated_snn_res.0.16 <- NULL 
total@meta.data$integrated_snn_res.0.17 <- NULL 
total@meta.data$integrated_snn_res.0.18 <- NULL 
total@meta.data$integrated_snn_res.0.19 <- NULL 
total@meta.data$integrated_snn_res.0.2 <- NULL 
total@meta.data$integrated_snn_res.0.21 <- NULL
total@meta.data$integrated_snn_res.0.22 <- NULL 
total@meta.data$integrated_snn_res.0.23 <- NULL 
total@meta.data$integrated_snn_res.0.24 <- NULL 
total@meta.data$integrated_snn_res.0.25 <- NULL 
total@meta.data$integrated_snn_res.0.26 <- NULL 
total@meta.data$integrated_snn_res.0.27 <- NULL
total@meta.data$integrated_snn_res.0.28 <- NULL 
total@meta.data$integrated_snn_res.0.29 <- NULL 
total@meta.data$integrated_snn_res.0.3 <- NULL 

###### Subset seurat object based on cell type and save seperate seurat objects per cell type

# CAMs
cell_type <- "CAMs"
Output.dir.subset = paste0(Output.dir, "6. Subcluster ",cell_type,"/")
dir.create(Output.dir.subset)
cells <- subset(total, Cell_category == cell_type)
cells
saveRDS(object = cells, file = paste0(Output.dir.subset, "6.0 Subset ", cell_type, ".rds"))

# Microglia without CAMs
cell_type <- "Microglia"
Output.dir.subset = paste0(Output.dir, "6. Subcluster ",cell_type,"/")
dir.create(Output.dir.subset)
cells <- subset(total, Cell_category == cell_type)
cells
saveRDS(object = cells, file = paste0(Output.dir.subset, "6.0 Subset ", cell_type, ".rds"))

# Microglia with CAMs
## Assign cell categories for subclustering
unique(cells$Cell_category)
## Alter cell category annotations to group Microglia with CAMs
total$Cell_category[total$Cell_category == "Microglia" | total$Cell_category == "CAMs"] = "Microglia+CAMs"
cell_type <- "Microglia+CAMs"
Output.dir.subset = paste0(Output.dir, "6. Subcluster ",cell_type,"/")
dir.create(Output.dir.subset)
cells <- subset(total, Cell_category == cell_type)
cells
saveRDS(object = cells, file = paste0(Output.dir.subset, "6.0 Subset ", cell_type, ".rds"))

# Astrocytes
cell_type <- "Astrocytes"
Output.dir.subset = paste0(Output.dir, "6. Subcluster ",cell_type,"/")
dir.create(Output.dir.subset)
cells <- subset(total, Cell_category == cell_type)
cells
saveRDS(object = cells, file = paste0(Output.dir.subset, "6.0 Subset ", cell_type, ".rds"))

# Oligodendrocytes
cell_type <- "Oligodendrocytes"
Output.dir.subset = paste0(Output.dir, "6. Subcluster ",cell_type,"/")
dir.create(Output.dir.subset)
cells <- subset(total, Cell_category == cell_type)
cells
saveRDS(object = cells, file = paste0(Output.dir.subset, "6.0 Subset ", cell_type, ".rds"))

# Neurons
cell_type <- "Neurons"
Output.dir.subset = paste0(Output.dir, "6. Subcluster ",cell_type,"/")
dir.create(Output.dir.subset)
cells <- subset(total, Cell_category == cell_type)
cells
saveRDS(object = cells, file = paste0(Output.dir.subset, "6.0 Subset ", cell_type, ".rds"))

# Immune
cell_type <- "Immune"
Output.dir.subset = paste0(Output.dir, "6. Subcluster ",cell_type,"/")
dir.create(Output.dir.subset)
cells <- subset(total, Cell_category == cell_type)
cells
saveRDS(object = cells, file = paste0(Output.dir.subset, "6.0 Subset ", cell_type, ".rds"))

# Endothelial
total$Cell_category[total$Cell_category == "Endothelial" | total$Cell_category == "Fibroblasts" 
                    | total$Cell_category == "Pericytes_vSMCs"] = "Endothelial+stromal"
cell_type <- "Endothelial+stromal"
Output.dir.subset = paste0(Output.dir, "6. Subcluster ",cell_type,"/")
dir.create(Output.dir.subset)
cells <- subset(total, Cell_category == cell_type)
cells
saveRDS(object = cells, file = paste0(Output.dir.subset, "6.0 Subset ", cell_type, ".rds"))

# Proceed to THR 6.1 Cell type subset Integration and clustering Round 1