# Donors have been genotyped for APOE, and have tau and amyloid pathology scores from the same tissue blocks used for snRNAseq. CellRanger count martices
# are stored in the input folder, with a seperate folder for each donor. The metadata .csv file is in the input folder. The output folder is created. 
# 3. Integration of samples for analysis
# Author: Thomas Rust (THR) adapted from Mirjam Koster (MKO)
# Date:23/02/2024 (updated)


# Installed modified version of Seurat for larger matrices. It modifies the spam matrix
#remotes::install_github("zhanghao-njmu/seurat")
# load packages
install.packages("Matrix", type = "source")
install.packages("irlba", type = "source")
library(Seurat)
library(reticulate)
library(dplyr)
library(patchwork)
library(ggplot2)

# Set output directory and read in Seurat object list with doublets removed by Scrublet
Output.dir <- "Output/"

total.list <- readRDS(file = paste0(Output.dir, "2.2 listed_samples doublets Scrublet removed.rds"))

#############################################################################
### Integration of multiple samples into 1 dataset
#############################################################################

# Determine variable features across samples for integration
total.features <- SelectIntegrationFeatures(object.list = total.list)

# Make sure all integration features are scaled
for(i in 1:length(total.list)){
  total.list[[i]] <- ScaleData(object = total.list[[i]], features = total.features)
  total.list[[i]] <- RunPCA(object = total.list[[i]], features = total.features)  
}

# Integrate datasets with rPCA 
# Decide the number of dimensions (by the elbow plots, usually 30 PCs is enough)
n.dims = 1:30
total.anchors <- FindIntegrationAnchors(object.list = total.list, reduction = 'rpca', dims = n.dims)

# save anchors
saveRDS(total.anchors, file = paste0(Output.dir, "3. anchors_integrated_samples.rds"))

#total.anchors <- readRDS(file = paste0(Output.dir, "3. anchors_integrated_samples.rds"))

n.dims = 1:30
total <- IntegrateData(anchorset = total.anchors, dims = n.dims, normalization.method = 'LogNormalize')
dim(total) #2000 x 38348
#DefaultAssay(total) <- "RNA" or "integrated"

saveRDS(total, file = paste0(Output.dir, "3. integrated_samples.rds"))
rm(list=setdiff(ls(), c("total", "Output.dir")))
gc()

# SCALING THE INTEGRATED DATA
# Linear transformation of expression (mean across cells = 0, variance across cells = 1)
# Scale all of the genes, not just the integrated ones
total <- readRDS(file = paste0(Output.dir, "3. integrated_samples.rds"))
all.genes <- rownames(total[["RNA"]])
total <- ScaleData(total, features = all.genes, vars.to.regress = c("percent.mito", "nCount_RNA")) # regress out percent mito RNA and number of transcripts
saveRDS(total, file = paste0(Output.dir, "3.1 integrated_dim_samples_scaled.rds"))
## DIMENSIONALITY REDUCTION OF THE INTEGRATED SAMPLES



## NON-LINEAR DIMENSIONAL REDUCTION

total <- readRDS(file = paste0(Output.dir, "3.1 integrated_dim_samples_scaled.rds"))

pdf(file = paste0(Output.dir, "3.1 Dimensionality reduction of integrated samples Scrublet doublets removed_PC30.pdf"))

# PCA of integrated samples to determine dimensions
total <- RunPCA(total, features = VariableFeatures(object = total)) #this is default
DimPlot(total, reduction = "pca", pt.size = 0.2, group.by = "sample", raster = TRUE)
ElbowPlot(total, ndims = 50)

#Decide the number of dimensions (of the integrated samples)
n.dims = 1:30
# UMAP, t-SNE
total <- RunUMAP(object = total, dims = n.dims)#, umap.method = "umap-learn", metric = "correlation")
total <- RunTSNE(object = total, dims = n.dims)
DimPlot(total, reduction = "umap", group.by = "sample", raster = TRUE)
DimPlot(total, reduction = "tsne", group.by = "sample", raster = TRUE)

dev.off()

# Examine and visualize PCA results using elbow plot and change number of PCs if needed

# Save the integrated data
saveRDS(total, file = paste0(Output.dir, "3.1 integrated_dim_samples Scrublet doublets removed_PC30.rds"))
rm(list=setdiff(ls(), c("total", "Output.dir")))
gc()

# Proceed to THR 4.0 Cluster for clustering of the integrated nuclei