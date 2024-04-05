# Analysis of snRNA-seq data
# 6.3 Subclustering of cell subsets starting with Round 2
# Author: Thomas Rust (adapted from Mirjam Koster)
# Date updated: 26/10/2023 (updated)

# load packages
library(Seurat)
library(ggplot2)
library(Matrix)
library(reshape2)
library(future)
library(clustree)

# set cell type
cell_type <- "Microglia"
# Options: Astrocytes Immune Neuron Oligodendocytes Endothelial+stromal Microglia+CAMs

# set output directory and subdirectory per cell type
Output.dir <- paste0("Output/")
Output.dir.subset = paste0(Output.dir, "6. Subcluster ",cell_type,"/CCA/")

##################  Load the subsetted Round 2 or Round 3 data #####################################################
Round = "Round 2" # CHANGE ROUND
cells <- readRDS(file = paste0(Output.dir.subset, "6.2.1 Input object", Round," ", cell_type, ".rds"))

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
  cells_intt <- IntegrateData(anchorset = anchors, dims = 1:30, normalization.method = "LogNormalize",
                            k.weight = 50)
} else {
  print(paste("The minimum number of nuclei per sample is", min(number_of_cells), "so k.weight set to this value"))
  cells_intt <- IntegrateData(anchorset = anchors, dims = 1:30, normalization.method = "LogNormalize", 
                             k.weight = min(number_of_cells))
}
#k.weight = 50 # solved an error by having non-default k.weight = 50 - Parameter depends on the nr of nuclei and library complexity

saveRDS(cells_intt, file = paste0(Output.dir.subset, "6.2.2 integrated ", Round, " ", integration.method,".rds"))

gc()
DefaultAssay(cells_intt) <- 'integrated'

# check if there are no NAs, sanity check. 
print('Sanity check: number of NAs in integrated data:' )
sum(is.na(cells_intt[["integrated"]]@data))
rm(list=setdiff(ls(), c("cells_intt", "Output.dir.subset","cell_type")))

###################################################################################################################
###                              Regressing & dimensionality reduction & plotting                               ###
###################################################################################################################

rm(list=setdiff(ls(), c("cells_intt", "Output.dir.subset","cell_type")))

# Scaling with regressed variables
cells_intt <- ScaleData(object = cells_intt, vars.to.regress = c("percent.mito", "percent.ribo", "nCount_RNA"))
#dataset dependent, check yourself what is relevant to regress out
gc()

Round <- "Round 2" # change round
pdf(file = paste0(Output.dir.subset, "6.2.2 Dimensionality reduction of ", cell_type," ", Round,".pdf"))

# PCA of integrated samples to determine dimensions
n.dims = 50
cells_intt <- RunPCA(object = cells_intt, features = VariableFeatures(object = cells_intt), npcs = n.dims) #this is default
DimPlot(cells_intt, reduction = "pca", pt.size = 0.2, group.by = "sample", raster = F)
ElbowPlot(cells_intt, ndims = n.dims) + 
  scale_x_continuous(breaks = seq(0,n.dims,10))

# Decide the number of dimensions (of the integrated samples). I chose 20 based off the ElbowPlot. It plateues 
# before 20 but I increased the PC to account for more variance 
n.dims = 1:20

# UMAP
cells_intt <- RunUMAP(object = cells_intt, dims = n.dims)
DimPlot(cells_intt, reduction = "umap", group.by = "sample", raster = F)
FeaturePlot(cells_intt, features = c("percent.mito","percent.ribo","nFeature_RNA","nCount_RNA"), 
            min.cutoff = "q1", max.cutoff = "q99", raster = F,
            cols = alpha(c("grey90", "#009DFF", "#C929DA","#D40000", "darkred"), 0.1))

dev.off()

# Save final object
saveRDS(cells_intt, file = paste0(Output.dir.subset, "6.2.2 Final integrated ", cell_type," ", Round,".rds"))

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
    print(VlnPlot(dataset, features = i, pt.size = F, raster = F) -
            FeaturePlot(object = dataset, features = i, label = TRUE, raster = F, reduction = "umap",
                        cols = alpha(c("#009DFF", "#C929DA","#D40000"), 0.5))
    )
  }
  dev.off()
}

###############################################################################################################################
###                                                 Subclustering                                                          ###
###############################################################################################################################

# Find neighbours
cells_int <- FindNeighbors(cells_int, dims = n.dims)

# Try multiple clustering resolutions to decide the optimal number of clusters
res <- seq(from = 0.0, to = 2, by = 0.1)
multiple.res.output(cells_int)

saveRDS(cells_int, file = paste0(Output.dir.subset, "6.2 Multiple res clustering ", cell_type," ", Round, ".rds")) 

############################## Final clustering/resolution ################################################################

res = 0.5 # EDIT based on results for this round
cells_int <- FindClusters(cells_int, resolution = res)
final.clustering.pdfs(cells_int)
ncol(cells_int)
table(cells_int$seurat_clusters)
# save final resubclustered object for this round
saveRDS(cells_int, file = paste0(Output.dir.subset, "6.2 ", cell_type, " ", Round, " resubclustered.rds"))

# Name the clusters
cells_int <- microglia

ID <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12") # change number of clusters manually based on the clustering resolution you select
names(ID) <- levels(cells_int)
cells_int$Clusters_subset <- Idents(cells_int)
cells_int <- RenameIdents(cells_int, ID)

Round <- "Round 3"
pdf(file = paste0(Output.dir.subset, "6.2 ", cell_type, " ", Round, " map plots.pdf"))
DimPlot(object = cells_int, reduction = "umap", label = T, raster = FALSE) + NoLegend()
DimPlot(object = cells_int, reduction = "umap", group.by = 'apoe', raster = FALSE)
DimPlot(object = cells_int, reduction = "umap", group.by = 'snRNAseq_batch', raster = FALSE)
DimPlot(object = cells_int, reduction = "umap", group.by = 'diagnosis', raster = FALSE )
    DimPlot(object = cells_int, reduction = "umap", split.by = 'diagnosis', raster = FALSE )
DimPlot(object = cells_int, reduction = "umap", group.by = 'group', raster = FALSE )
   
dev.off()

cells_int@meta.data$group <- as.factor(cells_int@meta.data$group)
cells_int@meta.data$apoe <- as.factor(cells_int@meta.data$apoe)
str(cells_int@meta.data$group)

Round <- "Round 3"
pdf(file = paste0(Output.dir.subset, "6.2 ", cell_type, " ", Round, " batch plots.pdf"))
DimPlot(object = cells_int, reduction = "umap", group.by = 'snRNAseq_batch', raster = FALSE)
DimPlot(object = cells_int, reduction = "umap", split.by = 'snRNAseq_batch', raster = FALSE)
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

ggplot(data = cells_int@meta.data, mapping = aes(x = apoe, fill = Clusters_subset)) +
  geom_bar() +
  theme_minimal() + 
  xlab("APOE (by grouped cluster)") +
  ylab("Cluster fraction") +
  ggtitle("Amount cluster per APOE genotype") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()