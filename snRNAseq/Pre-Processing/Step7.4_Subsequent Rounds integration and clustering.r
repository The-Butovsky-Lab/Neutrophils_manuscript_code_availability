# Analysis of snRNA-seq data
# 6.2.2 Subclustering of cell subsets starting with Round 2
# Author: Thomas Rust (adapted from Mirjam Koster)
# Date updated: 26/10/2023 (updated)

#remove.packages("Seurat")
# Installed modified version of Seurat for larger matrices. It modifies the spam matrix
remotes::install_github("zhanghao-njmu/seurat")

# load packages
library(Seurat)
library(ggplot2)
library(Matrix)
library(reshape2)
library(future)
library(clustree)
sessionInfo()

#remove.packages("SeuratObject")
#remotes::install_version(package = 'SeuratObject', version = package_version('4.1.3'))
#library("SeuratObject")
#sessionInfo()

# set cell type
cell_type <- "Microglia+CAMs"
# Options: Astrocytes Immune Neuron Oligodendocytes Endothelial+stromal Microglia+CAMs

# set output directory and subdirectory per cell type
Output.dir <- paste0("Output/")
Output.dir.subset = paste0(Output.dir, "6. Subcluster ",cell_type,"/")

##################  Load the subsetted Round 2 or Round 3 data #####################################################
Round = "Round 3" # CHANGE ROUND
cells <- readRDS(file = paste0(Output.dir.subset, "6.2.1 Input object ", Round," ", cell_type, ".rds"))

####################################################################################################################
###                                     Re-integration of subset                                                 ###
####################################################################################################################

total.number.of.nuclei = ncol(cells)

# Prepare the object for re-integration
DefaultAssay(cells) <- 'RNA'
cells[['integrated']]<- NULL #remove integrated assay
cells <- SplitObject(cells, split.by = 'sample') #split samples back into separate Seurat objects

print(length(cells))
gc()
for (i in 1:length(cells)) {
    print(i)
    cells[[i]] <- NormalizeData(cells[[i]], verbose = FALSE)
    cells[[i]] <- FindVariableFeatures(cells[[i]], selection.method = "vst", nfeatures = 2000, verbose = TRUE)
}

# Find anchors for CCA or rPCA
if(total.number.of.nuclei >= 200000){  
  #perform rPCA if more than 200,000 nuclei
  integration.method = "rPCA"
  print(paste0("Type of integration: ", integration.method, " because there's ", total.number.of.nuclei, " nuclei"))
  
  # Perform PCA per sample (only necessary for rPCA)
  features <- SelectIntegrationFeatures(object.list = cells, nfeatures = 2000)
  for (i in 1:length(cells)) {
    print(i)
    cells[[i]] <- ScaleData(object = cells[[i]], features = features)
    cells[[i]] <- RunPCA(object = cells[[i]], features = features)
  }
  
  # Find integration anchors
  anchors <- FindIntegrationAnchors(object.list = cells, dims = 1:30, reduction = "rpca",
                                    normalization.method = "LogNormalize", anchor.features = 2000)
  saveRDS(object = anchors, file = paste0(Output.dir.subset, "6.2.2 Anchors ", Round, " ",integration.method,".rds"))
    
} else { 
  # otherwise perform CCA
  integration.method = "CCA"
  print(paste0("Type of integration: ", integration.method, " because there's ", total.number.of.nuclei, " nuclei"))
  
  # Find integration anchors
  anchors <- FindIntegrationAnchors(object.list = cells, dims = 1:30, reduction = "cca", 
                                    normalization.method = "LogNormalize", anchor.features = 2000)
  saveRDS(object = anchors, file = paste0(Output.dir.subset, "6.2.2 Anchors ", Round, " ", integration.method,".rds"))
  
}
gc()

# read in anchors if there is an issue at integration step

# Perform the integration
sample.number.of.nuclei <- vector()
for(i in 1:length(cells)){
  sample.number.of.nuclei <- c(sample.number.of.nuclei, ncol(cells[[i]]))
}

#Adjust k.weight to minimun number of cells if necessary, otherwise k.weight = default 100
number_of_cells <- vector()

for(i in 1:length(cells)){
  number_of_cells <- c(number_of_cells, ncol(cells[[i]]))
}

if (min(number_of_cells) > 50){
  print(paste("The minimum number of nuclei per sample is", min(number_of_cells), "so default k.weight"))
  cells_int <- IntegrateData(anchorset = anchors, dims = 1:30, normalization.method = "LogNormalize",
                            k.weight = 50)
} else {
  print(paste("The minimum number of nuclei per sample is", min(number_of_cells), "so k.weight set to this value"))
  cells_int <- IntegrateData(anchorset = anchors, dims = 1:30, normalization.method = "LogNormalize", 
                             k.weight = min(number_of_cells))
}
#k.weight = 50 # solved an error by having non-default k.weight = 50 - Parameter depends on the nr of nuclei and library complexity

saveRDS(cells_int, file = paste0(Output.dir.subset, "6.2.2 integrated ", Round, " ", integration.method,".rds"))

gc()
DefaultAssay(cells_int) <- 'integrated'

# check if there are no NAs, sanity check. 
print('Sanity check: number of NAs in integrated data:' )
sum(is.na(cells_int[["integrated"]]@data))
#rm(list=setdiff(ls(), c("cells_int", "Output.dir.subset","cell_type")))

###################################################################################################################
###                              Regressing & dimensionality reduction & plotting                               ###
###################################################################################################################

# Scaling with regressed variables
cells_int <- ScaleData(object = cells_int, vars.to.regress = c("percent.mito", "percent.ribo", "nCount_RNA"))
#dataset dependent, check yourself what is relevant to regress out
gc()

Round <- "Round 3" # change round
pdf(file = paste0(Output.dir.subset, "6.2.2 Dimensionality reduction of ", cell_type," ", Round,".pdf"))

# PCA of integrated samples to determine dimensions
n.dims = 50
cells_int <- RunPCA(object = cells_int, features = VariableFeatures(object = cells_int), npcs = n.dims) #this is default
DimPlot(cells_int, reduction = "pca", pt.size = 0.2, group.by = "sample", raster = F)
ElbowPlot(cells_int, ndims = n.dims) + 
  scale_x_continuous(breaks = seq(0,n.dims,10))

# Decide the number of dimensions (of the integrated samples). I chose 20 based off the ElbowPlot. It plateues 
# before 20 but I increased the PC to account for more variance 
n.dims = 1:20

# UMAP
cells_int <- RunUMAP(object = cells_int, dims = n.dims)
DimPlot(cells_int, reduction = "umap", group.by = "sample", raster = F)
FeaturePlot(cells_int, features = c("percent.mito","percent.ribo","nFeature_RNA","nCount_RNA"), 
            min.cutoff = "q1", max.cutoff = "q99", raster = F,
            cols = alpha(c("grey90", "#009DFF", "#C929DA","#D40000", "darkred"), 0.1))

dev.off()

# Save final object
saveRDS(cells_int, file = paste0(Output.dir.subset, "6.2.2 Final integrated ", cell_type," ", Round,".rds"))

##########################################################################################################
###                                         Functions                                                  ###    
##########################################################################################################

# Known markers genes
known.markers0 <- c("P2RY12","CSF1R","C3","CX3CR1", #micro
                    "MRC1", #CAM
                    "GFAP","AQP4","SLC1A2","CABLES1", #astro
                    "PLP1","MOBP","MOG", #oligo, also: MBP
                    "PDGFRA","DSCAM","BCAN", #OPC
                    "RBFOX3","SLIT2", #neuron, also: "RELN"
                    "SLC17A7", #exci neuron
                    "GAD2","GRIP1", #inh neuron
                    "PTPRC", #immune
                    "IL7R","SKAP1","THEMIS", #T cells
                    "MS4A1","IGKC", #B cells
                    "NKG7","KLRD1",
                    "CLDN5","APOLD1","PECAM1", #endothelial, also CDH5
                    "PDGFRB","DCN", #stromal = pericyte + fibroblast
                    "DCDC1","CFAP299") #NK
# Some of these markers are from https://doi.org/10.1038/s41598-018-27293-5

known.markers <- paste("rna_", known.markers0, sep = "")

n.dims = 1:20

cols.continuous = alpha(c("grey90", "#009DFF", "#C929DA","#D40000", "darkred"), 0.3)

# Files for multiple resolution clustering
multiple.res.output <- function(dataset){
  
  res.clusters <- data.frame(res)
  
  pdf(paste0(Output.dir.subset, "6.2.2 ", cell_type, " ", Round, ".2 Clustering - Multiple res.pdf"), width = 16, height = 9)

  for(j in res){
    print(j)
    dataset <- FindClusters(dataset, resolution = j)
    res.clusters$Clustern[res.clusters$res == j] = nlevels(dataset)
    
    print(
      DimPlot(dataset, reduction = "umap", label = T, raster = F, shuffle = T) + ggtitle(paste0("UMAP, resolution: ",j)) -
            DimPlot(dataset, reduction = "umap", group.by = "sample", shuffle = T, raster = F)
    )
  }

  print("Vizualize cluster number and save")
  
  print(
    clustree(dataset)
  )
  
  print(
    ggplot(data = res.clusters, mapping = aes(x = res, y = Clustern)) +
      geom_point() +
      scale_x_continuous(breaks = res) +
      scale_y_continuous(breaks = seq(0, max(res.clusters$Clustern)+1, 1)) +
      theme_classic()
  )
  dev.off()

  write.table(x = res.clusters, file = paste0(Output.dir.subset, "6.2.2 ", cell_type, " ", Round, ".2 Clustering - Number of clusters.csv"), sep = ";")

}

# File for final clustering resolution
final.clustering.pdfs <- function(dataset) {
  
  pdf(file = paste0(Output.dir.subset, "6.2.2 ", cell_type, " ", Round, ".3 Clustering - Final res ", res,".pdf"), 
      width = 16, height = 9)
  print(
    DimPlot(dataset, reduction = "pca", label = TRUE, shuffle = T, raster = F) -
      DimPlot(dataset, reduction = 'pca', group.by = 'sample', shuffle = T, raster = F)
  )
  print(
    DimPlot(dataset, reduction = "umap", label = TRUE, shuffle = T, raster = F) -
      DimPlot(dataset, reduction = 'umap', group.by = 'sample', shuffle = T, raster = F)
  )
  
  print(
    VlnPlot(object = dataset, features = c("nCount_RNA","nFeature_RNA","percent.mito","percent.ribo"), 
            raster = F, pt.size = F, ncol = 3)
  )
  print(
    FeaturePlot(object = dataset, features = c("nCount_RNA","nFeature_RNA","percent.mito","percent.ribo"), 
                min.cutoff = 'q1', max.cutoff = 'q99', raster = F, ncol = 3, reduction = "umap",
                cols = cols.continuous)
  )
  
  print("Viualize known marker genes")
  print(
    DotPlot(object = dataset, assay = "RNA", features = known.markers0, cols = c("blue", "red"), dot.scale = 12) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  )
  
  for(i in known.markers0){
    print(VlnPlot(dataset, features = i, pt.size = F, raster = F, assay = "RNA") -
            FeaturePlot(object = dataset, features = i, label = TRUE, raster = F, reduction = "umap",
                        cols = alpha(c("#009DFF", "#C929DA","#D40000"), 0.5))
    )
  }
  dev.off()
}

###############################################################################################################################
###                                                 Subclustering                                                          ###
###############################################################################################################################

Round = "Round 3" # CHANGE ROUNDRound = "Round 2" # CHANGE ROUND

cells_int <- readRDS(file = paste0(Output.dir.subset, "6.2.2 Final integrated ", cell_type," ", Round,".rds"))

# Find neighbours
cells_int <- FindNeighbors(cells_int, dims = n.dims)

# Try multiple clustering resolutions to decide the optimal number of clusters
res <- seq(from = 0.0, to = 2, by = 0.1)
multiple.res.output(cells_int)

saveRDS(cells_int, file = paste0(Output.dir.subset, "6.2 Multiple res clustering ", cell_type," ", Round, ".rds")) 

############################## Final clustering/resolution ################################################################
Round <- "Round 3"
cells_int <- readRDS(file = paste0(Output.dir.subset, "6.2 Multiple res clustering ", cell_type," ", Round, ".rds")) 

# Microglia Round 2 - 0.6
# Microglia Round 3 - 0.6 
# Microglia+CAMs Round 2 - 0.6
# Microglia+CAMs Round 3 - 0.5

res = 0.5 # EDIT based on results for this round
cells_int <- FindClusters(cells_int, resolution = res)
final.clustering.pdfs(cells_int)
ncol(cells_int)
table(cells_int$seurat_clusters)
# save final resubclustered object for this round
saveRDS(cells_int, file = paste0(Output.dir.subset, "6.2 ", cell_type, " ", Round, " resubclustered.rds"))

# Name the clusters
# Microglia Round 2 - 14 clusters
# Microglia Round 3 - 13 clusters
# Microglia + CAMs Round 2 - 15 clusters 
# Microglia + CAMs Round 3 - 16 clusters 

ID <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15") # change number depdning on cell type object and round - resolution chosen of clusters manually based on the clustering resolution you select
names(ID) <- levels(cells_int)
cells_int$Clusters_subset <- Idents(cells_int)
cells_int <- RenameIdents(cells_int, ID)

Round <- "Round 3"
pdf(file = paste0(Output.dir.subset, "6.2 ", cell_type, " ", Round, " map plots.pdf"))
DimPlot(object = cells_int, reduction = "umap", label = T, raster = FALSE) + NoLegend()
DimPlot(object = cells_int, reduction = "umap", group.by = 'apoe', raster = FALSE)
DimPlot(object = cells_int, reduction = "umap", group.by = 'snRN.Aseq_batch', raster = FALSE)
DimPlot(object = cells_int, reduction = "umap", group.by = 'diagnosis', raster = FALSE )
    DimPlot(object = cells_int, reduction = "umap", split.by = 'diagnosis', raster = FALSE )
DimPlot(object = cells_int, reduction = "umap", group.by = 'group', raster = FALSE )
   
dev.off()

cells_int@meta.data$group <- as.factor(cells_int@meta.data$group)
cells_int@meta.data$apoe <- as.factor(cells_int@meta.data$apoe)
str(cells_int@meta.data$group)

pdf(file = paste0(Output.dir.subset, "6.2 ", cell_type, " ", Round, " batch plots.pdf"))
DimPlot(object = cells_int, reduction = "umap", group.by = 'snRN.Aseq_batch', raster = FALSE)
dev.off()

pdf(file = paste0(Output.dir.subset, "6.2.2 ", cell_type, " ", Round, " bar plots.pdf"))

ggplot(data = cells_int@meta.data, mapping = aes(x = Clusters_subset, fill = apoe)) +
  geom_bar(position = "fill") +
  theme_minimal() + 
  xlab("cluster") +
  ylab("Cell fraction") +
  ggtitle("Percentage of APOE genotypes per cluster") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = cells_int@meta.data, mapping = aes(x = Clusters_subset, fill = sample)) +
  geom_bar(position = "fill") +
  theme_minimal() + 
  xlab("cluster") +
  ylab("sample fraction") +
  ggtitle("Percentage of samples per cluster") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = cells_int@meta.data, mapping = aes(x = sample, fill = Clusters_subset)) +
  geom_bar(position = "fill") +
  theme_minimal() + 
  xlab("cluster") +
  ylab("sample fraction") +
  ggtitle("Percentage of cluster per sample") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = cells_int@meta.data, mapping = aes(x = Clusters_subset, fill = apoe)) +
  geom_bar() +
  theme_minimal() + 
  xlab("cluster") +
  ylab("Cell fraction") +
  ggtitle("Percentage of APOE genoptypes of total") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = cells_int@meta.data, mapping = aes(x = apoe, fill = Clusters_subset)) +
  geom_bar(position = "fill") +
  theme_minimal() + 
  xlab("APOE (by grouped cluster)") +
  ylab("Cluster fraction") +
  ggtitle("Amount cluster per apoe genotype") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = cells_int@meta.data, mapping = aes(x = apoe, fill = Clusters_subset)) +
  geom_bar() +
  theme_minimal() + 
  xlab("APOE (by grouped cluster)") +
  ylab("Cluster fraction") +
  ggtitle("Amount cluster per APOE genotype") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = cells_int@meta.data, mapping = aes(x = Clusters_subset, fill = diagnosis)) +
  geom_bar(position = "fill") +
  theme_minimal() + 
  xlab("cluster") +
  ylab("Cell fraction") +
  ggtitle("Percentage of AD/CTRL per cluster") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = cells_int@meta.data, mapping = aes(x = Clusters_subset, fill = sex)) +
  geom_bar(position = "fill") +
  theme_minimal() + 
  xlab("cluster") +
  ylab("Cell fraction") +
  ggtitle("Percentage of female/male per cluster") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = cells_int@meta.data, mapping = aes(x = Clusters_subset, fill = snRN.Aseq_batch)) +
  geom_bar(position = "fill") +
  theme_minimal() + 
  xlab("cluster") +
  ylab("Cell fraction") +
  ggtitle("Percentage of sequencing batch per cluster") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = cells_int@meta.data, mapping = aes(x = apoe, fill = Clusters_subset)) +
  geom_bar() +
  theme_minimal() +
  xlab("APOE (by grouped cluster)") +
  ylab("Cluster fraction") +
  ggtitle("Amount cluster per APOE genotype") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()

#### Microglia Round 3 

# Subset out smallest cluster for microglia and 2 smallest for microglia+CAMs and save Final objects for analysis

# switch to Round 3 for new data
Round = "Final"
cells_int
cells <- subset(x = cells_int, subset = seurat_clusters %in% c(0:13)) # change depedning on cell type
cells_int <- cells
ncol(cells_int)

# Name the clusters
ID <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13") # change number of clusters manually based on the clustering resolution you select
names(ID) <- levels(cells_int)
cells_int$Clusters_subset <- Idents(cells_int)
cells_int <- RenameIdents(cells_int, ID)

Round <- "Final"
pdf(file = paste0(Output.dir.subset, "6.2 ", cell_type, " ", Round, " map plots.pdf"))
DimPlot(object = cells_int, reduction = "umap", label = T, raster = FALSE) + NoLegend()
DimPlot(object = cells_int, reduction = "umap", group.by = 'apoe', raster = FALSE)
DimPlot(object = cells_int, reduction = "umap", group.by = 'snRN.Aseq_batch', raster = FALSE)
DimPlot(object = cells_int, reduction = "umap", group.by = 'diagnosis', raster = FALSE )
    DimPlot(object = cells_int, reduction = "umap", split.by = 'diagnosis', raster = FALSE )
DimPlot(object = cells_int, reduction = "umap", group.by = 'group', raster = FALSE )
   
dev.off()

cells_int@meta.data$group <- as.factor(cells_int@meta.data$group)
cells_int@meta.data$apoe <- as.factor(cells_int@meta.data$apoe)
str(cells_int@meta.data$group)

Round <- "Final"
pdf(file = paste0(Output.dir.subset, "6.2 ", cell_type, " ", Round, " batch plots.pdf"))
DimPlot(object = cells_int, reduction = "umap", group.by = 'snRN.Aseq_batch', raster = FALSE)
dev.off()

pdf(file = paste0(Output.dir.subset, "6.2.2 ", cell_type, " ", Round, " bar plots.pdf"))

ggplot(data = cells_int@meta.data, mapping = aes(x = Clusters_subset, fill = apoe)) +
  geom_bar(position = "fill") +
  theme_minimal() + 
  xlab("cluster") +
  ylab("Cell fraction") +
  ggtitle("Percentage of APOE genotypes per cluster") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = cells_int@meta.data, mapping = aes(x = Clusters_subset, fill = sample)) +
  geom_bar(position = "fill") +
  theme_minimal() + 
  xlab("cluster") +
  ylab("sample fraction") +
  ggtitle("Percentage of samples per cluster") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = cells_int@meta.data, mapping = aes(x = sample, fill = Clusters_subset)) +
  geom_bar(position = "fill") +
  theme_minimal() + 
  xlab("cluster") +
  ylab("sample fraction") +
  ggtitle("Percentage of cluster per sample") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = cells_int@meta.data, mapping = aes(x = Clusters_subset, fill = apoe)) +
  geom_bar() +
  theme_minimal() + 
  xlab("cluster") +
  ylab("Cell fraction") +
  ggtitle("Percentage of APOE genoptypes of total") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = cells_int@meta.data, mapping = aes(x = apoe, fill = Clusters_subset)) +
  geom_bar(position = "fill") +
  theme_minimal() + 
  xlab("APOE (by grouped cluster)") +
  ylab("Cluster fraction") +
  ggtitle("Amount cluster per apoe genotype") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = cells_int@meta.data, mapping = aes(x = apoe, fill = Clusters_subset)) +
  geom_bar() +
  theme_minimal() + 
  xlab("APOE (by grouped cluster)") +
  ylab("Cluster fraction") +
  ggtitle("Amount cluster per APOE genotype") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = cells_int@meta.data, mapping = aes(x = Clusters_subset, fill = diagnosis)) +
  geom_bar(position = "fill") +
  theme_minimal() + 
  xlab("cluster") +
  ylab("Cell fraction") +
  ggtitle("Percentage of AD/CTRL per cluster") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = cells_int@meta.data, mapping = aes(x = Clusters_subset, fill = sex)) +
  geom_bar(position = "fill") +
  theme_minimal() + 
  xlab("cluster") +
  ylab("Cell fraction") +
  ggtitle("Percentage of female/male per cluster") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = cells_int@meta.data, mapping = aes(x = Clusters_subset, fill = snRN.Aseq_batch)) +
  geom_bar(position = "fill") +
  theme_minimal() + 
  xlab("cluster") +
  ylab("Cell fraction") +
  ggtitle("Percentage of sequencing batch per cluster") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = cells_int@meta.data, mapping = aes(x = apoe, fill = Clusters_subset)) +
  geom_bar() +
  theme_minimal() +
  xlab("APOE (by grouped cluster)") +
  ylab("Cluster fraction") +
  ggtitle("Amount cluster per APOE genotype") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()

# save final resubclustered object for this round
Round <- "Final" # change round
saveRDS(cells_int, file = paste0(Output.dir.subset, "6.2 ", cell_type, " ", Round, " resubclustered.rds"))