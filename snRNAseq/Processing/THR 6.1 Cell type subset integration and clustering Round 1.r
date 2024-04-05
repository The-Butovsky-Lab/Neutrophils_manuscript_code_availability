# Analysis of snRNA-seq data
# 6.1 Integration of cell subsets with CCA 
# Author: Thomas Rust adapted from Mirjam Koster
# Date updated: 23/02/2024

# load packages
library(Seurat)
library(ggplot2)
library(Matrix)
library(reshape2)
library(future)

# Set cell type 
cell_type <- "Microglia+CAMs"
# Options: Astrocytes Immune Neuron Oligodendocytes Endothelial+stromal Microglia+CAMs Microglia

# set output directory and subdirectory
Output.dir <- paste0("Output/")
Output.dir.subset = paste0(Output.dir, "6. Subcluster ",cell_type,"/")

# Load the subsetted data
cells <- readRDS(paste0(file = Output.dir.subset, "6.0 Subset ", cell_type, ".rds"))

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
  #perform rPCA if more than 200,000 nuclei, all clusters are below 200k
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
  saveRDS(object = anchors, file = paste0(Output.dir.subset, "6.1.1 Anchors ", integration.method,".rds"))
    
} else { 
  # otherwise perform CCA
  integration.method = "CCA"
  print(paste0("Type of integration: ", integration.method, " because there's ", total.number.of.nuclei, " nuclei"))
  
  # Find integration anchors
  anchors <- FindIntegrationAnchors(object.list = cells, dims = 1:30, reduction = "cca", 
                                    normalization.method = "LogNormalize", anchor.features = 2000)
  saveRDS(object = anchors, file = paste0(Output.dir.subset, "6.1.1 Anchors ", integration.method,".rds"))
  
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
  cells_int <- IntegrateData(anchorset = anchors, dims = 1:30, normalization.method = "LogNormalize",
                            k.weight = 50)
} else {
  print(paste("The minimum number of nuclei per sample is", min(number_of_cells), "so k.weight set to this value"))
  cells_int <- IntegrateData(anchorset = anchors, dims = 1:30, normalization.method = "LogNormalize", 
                             k.weight = min(number_of_cells))
}
#k.weight = 50 # solved an error by having non-default k.weight = 50 - Parameter depends on the nr of nuclei and library complexity

saveRDS(cells_int, file = paste0(Output.dir.subset, "6.1.2 integrated_cca.rds"))

gc()
DefaultAssay(cells_int) <- 'integrated'

# check if there are no NAs, sanity check. 
print('Sanity check: number of NAs in integrated data:' )
sum(is.na(cells_int[["integrated"]]@data))
#rm(list=setdiff(ls(), c("cells_int", "Output.dir.subset","cell_type")))

##################################################################################################################
###                              Regressing & dimensionality reduction & plotting                              ###
###################################################################################################################

# Scaling with regressed variables
cells_int <- ScaleData(object = cells_int, vars.to.regress = c("percent.mito", "percent.ribo", "nCount_RNA")) #dataset dependent, check yourself what is relevant to regress out
gc()

pdf(file = paste0(Output.dir.subset, "6.1.3 Dimensionality reduction of ", cell_type,".pdf"))

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
saveRDS(cells_int, file = paste0(Output.dir.subset, "6.1.3 Final integrated ", cell_type,".rds"))

###################################################################################################################
###                                           Subclustering round 1                                             ###
###################################################################################################################

Cell_subset <- cells_int
Round <- "Round 1"

# Find neighbours
Cell_subset <- FindNeighbors(Cell_subset, dims = n.dims)

# Try multiple clustering resolutions to decide the optimal number of clusters
res <- seq(from = 0.0, to = 2, by = 0.1)
res.clusters <- data.frame(res)

pdf(paste0(Output.dir.subset, "6.1 ", Round,".1 Clustering - res 0 to 2 by 0.1.pdf"), 
    width = 16, height = 9)

for(j in res){
  print(j)
  Cell_subset <- FindClusters(Cell_subset, resolution = j)
  res.clusters$Clustern[res.clusters$res == j] = nlevels(Cell_subset)
  print(DimPlot(Cell_subset, reduction = "umap", label = T ,shuffle = T, raster = F) + 
          ggtitle(paste0("UMAP, resolution: ",j)) -
          DimPlot(Cell_subset, reduction = "umap", group.by = "sample", shuffle = T, raster = F))
}

dev.off()

pdf(paste0(Output.dir.subset, "6.1 ", Round,".2 Clustering - number of clusters.pdf"), width = 16, height = 9)

ggplot(data = res.clusters, mapping = aes(x = res, y = Clustern)) +
  geom_point() +
  scale_x_continuous(breaks = res) +
  scale_y_continuous(breaks = seq(0, max(res.clusters$Clustern)+1, 1)) +
  theme_classic()
dev.off()

# save Round 1 object at multiple res
saveRDS(Cell_subset, file = paste0(Output.dir.subset, "6.1.4 Mulitple res ", cell_type,".rds"))
write.table(x = res.clusters, file = paste0(Output.dir.subset, "6.1 ", Round,".3 Clustering - number of clusters.csv"), sep = ";")

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


################# Final clustering resolution - Round 1
# Function for final clustering resolution
cols.continuous = alpha(c("grey90", "#009DFF", "#C929DA","#D40000", "darkred"), 0.3)

final.clustering.pdfs <- function(dataset) {
  
  pdf(file = paste0(Output.dir.subset, "6.1 ", cell_type, " ", Round, ".3 Clustering - Final res ", res,".pdf"), 
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
    DotPlot(object = dataset, features = known.markers0, cols = c("blue", "red"), assay = "RNA",dot.scale = 12) +
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

############################################

#################### Final clustering/resolution ################################

# set clustering to round 1
Round <- "Round 1"

# read in cell type object with multiple clustering resolutions
cells <-readRDS(file = paste0(Output.dir.subset, "6.1.4 Mulitple res ", cell_type,".rds"))

# set final resolution:
# Microglia 0.5
#Microglia+CAMs 0.7

res = 0.7 # EDIT based on results
cells <- FindClusters(cells, resolution = res)
final.clustering.pdfs(cells)

ncol(cells)
table(cells$seurat_clusters)

########### Name the clusters############################
# Microglia - 18 clusters
# Microglia + CAMs - 20 clusters
ID <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19") # change number of clusters manually each time
names(ID) <- levels(cells)
cells$Clusters_subset <- Idents(cells)
cells <- RenameIdents(cells, ID)

# Create output files
pdf(file = paste0(Output.dir.subset, "6.1 ", cell_type," ", Round, " map plots.pdf"))
DimPlot(object = cells, reduction = "umap", label = T, raster = FALSE) + NoLegend()
DimPlot(object = cells, reduction = "umap", group.by = 'apoe', raster = FALSE)
DimPlot(object = cells, reduction = "umap", split.by = 'apoe', raster = FALSE )
DimPlot(object = cells, reduction = "umap", group.by = 'diagnosis', raster = FALSE )
    DimPlot(object = cells, reduction = "umap", split.by = 'diagnosis', raster = FALSE )
DimPlot(object = cells, reduction = "umap", group.by = 'group', raster = FALSE )
DimPlot(object = cells, reduction = "umap", group.by = 'snRN.Aseq_batch', raster = FALSE )
   
dev.off()

cells@meta.data$group <- as.factor(cells@meta.data$group)
cells@meta.data$apoe <- as.factor(cells@meta.data$apoe)
str(cells@meta.data$group)

# save plots in pdf format
pdf(file = paste0(Output.dir.subset, "6.1 ", cell_type," ", Round, " bar plots.pdf"))

ggplot(data = cells@meta.data, mapping = aes(x = Clusters_subset, fill = apoe)) +
  geom_bar(position = "fill") +
  theme_minimal() + 
  xlab("cluster") +
  ylab("Cell fraction") +
  ggtitle("Percentage of APOE genotypes per cluster") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = cells@meta.data, mapping = aes(x = Clusters_subset, fill = sample)) +
  geom_bar(position = "fill") +
  theme_minimal() + 
  xlab("cluster") +
  ylab("sample fraction") +
  ggtitle("Percentage of samples per cluster") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = cells@meta.data, mapping = aes(x = sample, fill = Clusters_subset)) +
  geom_bar(position = "fill") +
  theme_minimal() + 
  xlab("cluster") +
  ylab("sample fraction") +
  ggtitle("Percentage of cluster per sample") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = cells@meta.data, mapping = aes(x = Clusters_subset, fill = apoe)) +
  geom_bar() +
  theme_minimal() + 
  xlab("cluster") +
  ylab("Cell fraction") +
  ggtitle("Percentage of APOE genoptypes of total") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = cells@meta.data, mapping = aes(x = apoe, fill = Clusters_subset)) +
  geom_bar(position = "fill") +
  theme_minimal() + 
  xlab("APOE (by grouped cluster)") +
  ylab("Cluster fraction") +
  ggtitle("Amount cluster per apoe genotype") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = cells@meta.data, mapping = aes(x = apoe, fill = Clusters_subset)) +
  geom_bar() +
  theme_minimal() + 
  xlab("APOE (by grouped cluster)") +
  ylab("Cluster fraction") +
  ggtitle("Amount cluster per APOE genotype") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = cells@meta.data, mapping = aes(x = Clusters_subset, fill = diagnosis)) +
  geom_bar(position = "fill") +
  theme_minimal() + 
  xlab("cluster") +
  ylab("Cell fraction") +
  ggtitle("Percentage of AD/CTRL per cluster") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = cells@meta.data, mapping = aes(x = Clusters_subset, fill = sex)) +
  geom_bar(position = "fill") +
  theme_minimal() + 
  xlab("cluster") +
  ylab("Cell fraction") +
  ggtitle("Percentage of female/male per cluster") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = cells@meta.data, mapping = aes(x = Clusters_subset, fill = snRN.Aseq_batch)) +
  geom_bar(position = "fill") +
  theme_minimal() + 
  xlab("cluster") +
  ylab("Cell fraction") +
  ggtitle("Percentage of sequencing batch per cluster") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = cells@meta.data, mapping = aes(x = apoe, fill = Clusters_subset)) +
  geom_bar() +
  theme_minimal() + 
  xlab("APOE (by grouped cluster)") +
  ylab("Cluster fraction") +
  ggtitle("Amount cluster per APOE genotype") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()

# save the Round 1 objects for the reintegrated and reclustered objects for each cell type
saveRDS(cells, file = paste0(Output.dir.subset, "6.1 ", cell_type, " ", Round, " resubclustered.rds"))

# Proceed to THR 6.2 Subsequent Round of cell type integration and clustering