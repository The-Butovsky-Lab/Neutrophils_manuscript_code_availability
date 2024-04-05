# Donors have been genotyped for APOE, and have tau and amyloid pathology scores from the same tissue blocks used for snRNAseq. CellRanger count martices
# are stored in the input folder, with a seperate folder for each donor. The metadata .csv file is in the input folder. The output folder is created. 
# 4. Clustering of nuclei in the integrated dataset
# Author: Thomas Rust (THR) adapted from Mirjam Koster (MKO)
# Date: 26/10/2023 (updated)

# load packages
library(Seurat)
library(reticulate)
library(dplyr)
library(patchwork)
library(ggplot2)
library(clustree)

# set output directory and read in integrated dataset
Output.dir <- "Output/"

total.list <- readRDS(file = paste0(Output.dir, "3.1 integrated_dim_samples Scrublet doublets removed_PC30.rds"))

#############################################################################
### Clustering
#############################################################################

# set number of dimensions and Find Neighbours for clustering 
n.dims = 1:30
total.list <- FindNeighbors(total.list, dims = n.dims)

#To decide the optimal number of clusters, try different resolutions in a rnage from 0.0 to 0.3, in increments of 0.01
res <- seq(from = 0.0, to = 0.3, by = 0.01)
res.clusters <- data.frame(res)

# save plots for mulitple resolution clustering analysis
pdf(paste0(Output.dir, "4.1 Clustering - multiple resolution.pdf"), 
    width = 16, height = 9)

for(j in res){
  total.list <- FindClusters(total.list, resolution = j)
  res.clusters$Clustern[res.clusters$res == j] = nlevels(total.list)
  
  print(DimPlot(total.list, reduction = "pca", label = TRUE, raster = TRUE) + ggtitle(paste0("PCA, resolution: ",j)) - 
          DimPlot(total.list, reduction = "pca", group.by = "sample", raster = TRUE))
  print(DimPlot(total.list, reduction = "umap", label = TRUE, raster = TRUE) + ggtitle(paste0("UMAP, resolution: ",j)) -
          DimPlot(total.list, reduction = "umap", group.by = "sample", raster = TRUE))
  print(DimPlot(total.list, reduction = "tsne", label = TRUE, raster = TRUE) + ggtitle(paste0("TSNE, resolution: ",j)) - 
          DimPlot(total.list, reduction = "tsne", group.by = "sample", raster = TRUE))
}

dev.off()

# save plots for mulitple resolution clustering analysis
pdf(paste0(Output.dir, "4.1 Clustering - cluster number.pdf"), width = 16, height = 9)
print(clustree(total.list))
ggplot(data = res.clusters, mapping = aes(x = res, y = Clustern)) +
  geom_point() +
  scale_x_continuous(breaks = res) +
  scale_y_continuous(breaks = seq(0, max(res.clusters$Clustern)+1, 1)) +
  theme_classic()
dev.off()

# save number of clusters per resolution as .csv file for mulitple resolution clustering analysis
write.table(x = res.clusters, file = paste0(Output.dir, "4.1 Clustering - multiple res.csv"), sep = ";")

# save the seurat object with the nuclei assigned to a cluster at different cluster resolutions, of which we will choose one
saveRDS(total.list, file = paste0(Output.dir, "4.0 all_samples_clustered_multiple_res.rds"))

#### END OF SEMI-AUTOMATIC PART OF THE CODE
# Please do check if all outputs appear normal, particularly:
# 2. Quality control: check if 5% mito was a realistic cut-off
# 3. Normalised data: does it look nicely bell-curved?
# 5. PCA: is 30 PCs enough to cover most variability of the individual samples?
# 6. PCA: is 30 PCs enough to cover most variability of the integrated samples?
# 7. Resolutions: determine the number of cluster you would like to keep for the next step

# Check object
total.list

#Final clustering/resolution - choose optimal resolution for your dataset, normally where there is a plateua in the resolotion by cluster number plot
res = 0.25 # EDIT based on results
total.list <- FindClusters(total.list, resolution = res)

pdf(file = paste0(Output.dir, "4.2 Clustering - final res", res,".pdf"), 
    width = 16, height = 9)
DimPlot(total.list, reduction = "pca", label = TRUE, raster = TRUE) -
  DimPlot(total.list, reduction = 'pca', group.by = 'sample', raster = TRUE)
DimPlot(total.list, reduction = "umap", label = TRUE, raster = TRUE) -
  DimPlot(total.list, reduction = 'umap', group.by = 'sample', raster = TRUE)
DimPlot(total.list, reduction = "tsne", label = TRUE, raster = TRUE) -
  DimPlot(total.list, reduction = 'tsne', group.by = 'sample', raster = TRUE)

dev.off()


## IDENTIFYING THE CLUSTERS 
#Expression of known markers genes
known.markers0 <- c("P2RY12","CSF1R","C3","CX3CR1", #microglia
                    "MRC1", #CAMs
                    "GFAP","AQP4","SLC1A2","CABLES1", #astrocytes
                    "PLP1","MOBP","MOG", "MBP", #oligo, also: MBP
                    "PDGFRA","DSCAM","BCAN","VCAN", #OPC
                    "RBFOX3","SLIT2", #neuron, also: "RELN"
                    "SLC17A7", #exci neuron
                    "GAD2","GRIP1", #inh neuron
                    "PTPRC", #immune
                    "IL7R","SKAP1","THEMIS", #T cells
                    "MS4A1","IGKC", #B cells
                    "NKG7","KLRD1", # NK cells
                    "CLDN5","APOLD1","PECAM1", #endothelial, also CDH5
                    "PDGFRB","DCN", #stromal = pericyte + fibroblast, also PDGFRA
                    "CALD1", # smooth muscle cells
                    "DCDC1","CFAP299") #ependymal
# Some of these markers are from https://doi.org/10.1038/s41598-018-27293-5

known.markers <- paste("rna_", known.markers0, sep = "")

# save pdf of known marker expression
pdf(file = paste0(Output.dir, "4.3 Known marker genes_raster.pdf"),
    width = 16, height = 9)
DotPlot(object = total.list, assay = "RNA", features = known.markers0, cols = c("blue", "red"), dot.scale = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
for(i in known.markers){
  print(VlnPlot(total.list, features = i, pt.size = 0) -
          FeaturePlot(object = total.list, features = i, reduction = "umap", label = TRUE, raster = TRUE,
                      cols = alpha(c("#009DFF", "#C929DA","#D40000"), 0.5))
  )
}
dev.off()

## NAME THE CLUSTERS
# EDIT based on results - broad neuron clusters 
ID <- c("Astrocytes1", #0
        "Microglia1", #1
        "Microglia2", #2
        "Astrocytes2", #3
        "Astrocytes3", #4
        "Oligodendrocytes1", #5
        "Endothelial1", #6
        "Pericytes_vSMCs", #7
        "Neurons1", #8
        "Immune1", #9
        "Fibroblasts1", #10
        "Neurons2", #11
        "CAMs", # 12
        "Neurons3", #13
        "Immune2", #14
        "Neurons4", #15
        "Neurons5", #16
        "Neurons6", #17
        "Astrocytes4") #18
names(ID) <- levels(total.list)

# rename idents to cell types
total.list$Clusters_n <- Idents(total.list)
total.list <- RenameIdents(total.list, ID)
total.list$Cell_types <- Idents(total.list)

# save the seurat object with nuclei martked with their individual cell type clusters
saveRDS(total.list, file = paste0(Output.dir, "4.2 all_samples.rds"))
rm(list=setdiff(ls(), c("total.list", "Output.dir")))

#Final visualization for the seurat object with nuclei martked with their individual cell type clusters
total <- total.list

pdf(file = paste0(Output.dir, "4.5.Final visualisations.pdf"))
DimPlot(object = total, reduction = "umap", group.by = "Cell_types", label = TRUE, raster = TRUE) + NoLegend()
DimPlot(object = total, reduction = "umap", group.by = "sample", raster = TRUE)
DimPlot(object = total, reduction = "tsne", group.by = "Cell_types", label = TRUE, raster = TRUE) + NoLegend() 
DimPlot(object = total, reduction = "tsne", group.by = "sample", raster = TRUE)
DimPlot(object = total, reduction = "pca", group.by = "Cell_types", label = TRUE, repel = TRUE, raster = TRUE) + NoLegend()
DimPlot(object = total, reduction = "pca", group.by = "sample", raster = TRUE)
VlnPlot(object = total, features = c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo"), 
        pt.size = 0, ncol = 2, assay = 'RNA') & theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
VlnPlot(object = total, features = c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo"), 
        pt.size = 0, ncol = 2, assay = 'RNA', group.by = "sample") & theme(axis.text.x = element_text(size = 4, angle = 90, vjust = 0.5, hjust=1))

ggplot(data = total@meta.data, mapping = aes(x = sample, fill = Cell_types)) +
  geom_bar(position = "fill") +
  theme_minimal() + 
  xlab("Sample") +
  ylab("Cell fraction") +
  ggtitle("Percentage of cell types of total") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = total@meta.data, mapping = aes(x = sample, fill = Cell_types)) +
  geom_bar() +
  theme_minimal() + 
  xlab("Sample") +
  ylab("Cell fraction") +
  ggtitle("Percentage of cell types of total") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = total@meta.data, mapping = aes(x = Cell_types, fill = sample)) +
  geom_bar(position = "fill") +
  theme_minimal() + 
  xlab("Cell type (by cluster)") +
  ylab("Sample fraction") +
  ggtitle("Percentage sample per cluster") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

## OPTIONAL: COMBINE CLUSTERS PER CELL TYPE

total <- readRDS(file = paste0(Output.dir, "4.2 all_samples.rds"))
total@meta.data

# If there's multiple clusters for 1 or more cell types it is possible to 
# combine those to get a better idea of cell type distribution
total$Cell_category <- gsub(pattern = "[0-9, ?]", replacement = "", x = total$Cell_types)

pdf(file = paste0(Output.dir, "4.6.Final visualisations - combined clusters.pdf"))
DimPlot(object = total, reduction = "umap", group.by = "Cell_category", label = T, raster = TRUE) + NoLegend()
DimPlot(object = total, reduction = "umap", group.by = "sample", raster = TRUE) + NoLegend()
DimPlot(object = total, reduction = "umap", group.by = "diagnosis", raster = TRUE) + NoLegend()
DimPlot(object = total, reduction = "umap", group.by = "apoe", raster = TRUE) + NoLegend()

ggplot(data = total@meta.data, mapping = aes(x = sample, fill = Cell_category)) +
  geom_bar(position = "fill") +
  theme_minimal() + 
  xlab("Sample") +
  ylab("Cell fraction") +
  ggtitle("Percentage of cell types of total") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = total@meta.data, mapping = aes(x = sample, fill = Cell_category)) +
  geom_bar() +
  theme_minimal() + 
  xlab("Sample") +
  ylab("Cell count") +
  ggtitle("Amount of cell types of total") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = total@meta.data, mapping = aes(x = Cell_category, fill = sample)) +
  geom_bar(position = "fill") +
  theme_minimal() + 
  xlab("Cell type (by grouped cluster)") +
  ylab("Sample fraction") +
  ggtitle("Percentage sample per cluster") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = total@meta.data, mapping = aes(x = Cell_category, fill = sample)) +
  geom_bar() +
  theme_minimal() + 
  xlab("Cell type (by grouped cluster)") +
  ylab("Sample fraction") +
  ggtitle("Percentage sample per cluster") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()

# Calculate cell proportions
cell.num <- table(total@meta.data$Cell_category)
cell.perc <-  round((cell.num/sum(cell.num))*100, digits = 1)
cell.num
cell.perc
Cell_category_fraction_table <- table(total@meta.data$Cell_category, total@meta.data$group)
Cell_category_fraction_table

# Add cell number per cluster to cluster labels
#total <- readRDS(file = paste0(Output.dir, "4.2 all_samples_cell_category.rds"))

ClusterLabels = paste("Cluster:", names(cell.perc), paste0("(", cell.perc,"%)"))

ClusterLabels = paste(names(cell.perc), paste0(": ", cell.perc, "%")) 

ClusterLabelsOrdered <- c("Astrocytes : 46.4%", "Microglia : 34.5%", "Oligodendrocytes 5.9:%", "Neurons 4.2: %", 
                          "Endothelial : 3.0%", "Pericytes_vSMCs : 2.5%", "Immune : 1.5%", "Fibroblasts : 1.0%", "CAMs : 0.9%")

# Order legend labels in plot in the same order as 'ClusterLabels'

ClusterPerc = names(cell.perc)
cols.continuous = alpha(c("grey90", "#009DFF", "#C929DA","#D40000", "darkred"), 0.3)
pdf(file = paste0(Output.dir, "4.6.1 Final visualisations - combined clusters.pdf"))
DimPlot(object = total, reduction = "umap", group.by = "Cell_category", label = T, raster = TRUE) +
  scale_colour_discrete(breaks = ClusterPerc, 
                        labels = ClusterLabelsOrdered) + labs(x = "UMAP 1", y = "UMAP 2")
VlnPlot(object = total, features = c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo"), 
        pt.size = 0, ncol = 2, assay = 'RNA', group.by = "Cell_category") & theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
FeaturePlot(total, features = c("percent.mito","percent.ribo","nFeature_RNA","nCount_RNA"), 
              min.cutoff = "q1", max.cutoff = "q99", raster = F, pt.size = 1/10^6, order = T)
dev.off()

write.csv2(x = table(total@meta.data$Cell_category), file = paste0(Output.dir,"4.7 Cell type numbers_cell_category.csv"))
Cell_category_fraction_table <- table(total@meta.data$Cell_types, total@meta.data$sample)
Cell_category_fraction_table
write.csv2(x = table(total@meta.data$Cell_category, total@meta.data$sample), 
           file = paste0(Output.dir,"4.8 Cell type numbers by sample.csv"))
Cell_category_fraction_table <- table(total@meta.data$Cell_category, total@meta.data$group)
Cell_category_fraction_table
write.csv2(x = table(total@meta.data$Cell_category, total@meta.data$group), 
           file = paste0(Output.dir,"4.9 Cell type numbers by group.csv"))
Cell_category_fraction_table <- table(total@meta.data$Cell_category, total@meta.data$diagnosis)
Cell_category_fraction_table
write.csv2(x = table(total@meta.data$Cell_category, total@meta.data$diagnosis), 
           file = paste0(Output.dir,"4.9 Cell type numbers by diagnosis.csv"))

# SAVE THE SEURAT OBJECT WITH NUCLEI GROUPED IN TO CELL TYPES - use this for analysis
saveRDS(total, file = paste0(Output.dir, "4.2 all_samples_cell_category.rds"))

# Proceed to THR 5.0 Update object metadata