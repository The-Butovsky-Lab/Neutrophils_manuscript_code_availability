# Analysis of snRNA-seq data
# 6.2.1 Subsetting of sub clusters for reclustering, starting with Round 2
# After Round 1 of sublcutering, remove clusters that are either doublets, predominantly from one or town donors, or heavily infulenced by
# mitochondrial or ribosomal RNA. Remove those clusters from the subset using this script and save a new input object. Run the input object through 'Round 2' 
# reintegration and clustering. Keep doing this until you have clean clusters for each cell type for downstream analysis. 
# Author: Thomas Rust
# Date updated: 23/02/2024 (updated)

# load packages
library(Seurat)
library(ggplot2)
library(Matrix)
library(reshape2)
library(future)
library(clustree)

# set cell type
cell_type <- "Microglia+CAMs"
# Options: Astrocytes Immune Neuron Oligodendocytes Endothelial+stromal Microglia+CAMs

# set output directory and subdirectory per cell type
Output.dir <- paste0("Output/")
Output.dir.subset = paste0(Output.dir, "6. Subcluster ",cell_type,"/")

sessionInfo()

# Load the subclustered Round 1 data 
Round = "Round 1"
cells <- readRDS(file = paste0(Output.dir.subset, "6.1 ", cell_type, " ", Round, " resubclustered.rds"))

###############################################################################################################################
###                                                Microglia Subclustering round 2                                                   ###
###############################################################################################################################
# remove 13 and 17 due to abundance of a few samples and possible batch effect, mito content
# remove 12 and 16 as they appear to be microglia-astrocyte doublets
# remove cluster 6 and 10 due to ribo content
# switch to Round 2 for new data
Round = "Round 2"

cells2 <- subset(x = cells, subset = seurat_clusters %in% c(0:5,7:9,11,14,15))
ncol(cells2)

saveRDS(cells2, file = paste0(Output.dir.subset, "6.2.1 Input object ", Round," ", cell_type, ".rds"))

# Load the subclustered Round 2 data 
Round = "Round 2"
cells <- readRDS(file = paste0(Output.dir.subset, "6.2 ", cell_type, " ", Round, " resubclustered.rds"))

###############################################################################################################################
###                                                 Microglia Subclustering round 3                                                   ###
###############################################################################################################################
# remove cluster 11 due to high ribo content and low counts
# switch to Round 3 for new data
Round = "Round 3"

cells2 <- subset(x = cells, subset = seurat_clusters %in% c(0:10,12,13))
ncol(cells2)

saveRDS(cells2, file = paste0(Output.dir.subset, "6.2.1 Input object ", Round," ", cell_type, ".rds"))

###############################################################################################################################
###                                                Microglia+CAMs Subclustering round 2                                                   ###
###############################################################################################################################
# remove 6 and 13 due to ribo content
# remove 14 as they appear to be microglia-astrocyte doublets
# remove cluster 18 and 19 due to mito content
# switch to Round 2 for new data
Round = "Round 2"

cells2 <- subset(x = cells, subset = seurat_clusters %in% c(0:5,7:12,15:17))
ncol(cells2)

saveRDS(cells2, file = paste0(Output.dir.subset, "6.2.1 Input object ", Round," ", cell_type, ".rds"))

# Load the subclustered Round 2 data 
Round = "Round 2"
cells <- readRDS(file = paste0(Output.dir.subset, "6.2 ", cell_type, " ", Round, " resubclustered.rds"))

###############################################################################################################################
###                                                 Microglia+CAMs Subclustering round 3                                                   ###
###############################################################################################################################
# remove cluster 11 due to high ribo content and low counts
# switch to Round 3 for new data
Round = "Round 3"

cells2 <- subset(x = cells, subset = seurat_clusters %in% c(0:12))
ncol(cells2)

saveRDS(cells2, file = paste0(Output.dir.subset, "6.2.1 Input object ", Round," ", cell_type, ".rds"))