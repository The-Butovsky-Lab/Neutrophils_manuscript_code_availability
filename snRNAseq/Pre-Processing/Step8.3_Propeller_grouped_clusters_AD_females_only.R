# 7.3. Propeller for calculating cell type proportions
# Test metadata variables for effects on cell type proportions
# Test these parameters in grouped clusters for AD females only
# Date: 19/02/2024

# load packages
library(Seurat)
library(ggplot2)
library(data.table)
library(readxl)
library(rstatix)
library(dplyr)
library(ggpubr)
library(speckle)
library(statmod)
library(limma)
library(gt)

# load in microglia seurat object
cell_type = "Microglia"
Round <- "Round 3"
setwd("D:/Dropbox/snRNAseq")
Output.dir.subset = paste0("D:/Dropbox/snRNAseq/Neta_APOE33_34_analysis/results/Microglia/Cluster_proportions/Grouped_clusters/")
Microglia_subset <- readRDS(file = paste0("D:/Dropbox/snRNAseq/Neta_APOE33_34_analysis/Object/6.2 MicrogliaFinalresubclustered_updated-ClusterNames.rds"))

# Name the clusters
ID <- c("1", "3", "8", "10", "4", "9", "2", "7", "5", "6", "11", "12") #### change number of clusters
names(ID) <- levels(Microglia_subset)
levels(Microglia_subset)
Microglia_subset$reordered_mappings <- Idents(Microglia_subset)
Microglia_subset <- RenameIdents(Microglia_subset, ID)
levels(Microglia_subset)

## Assign cell categories for subclustering
# Group homeostatic clusters 1 and 2 group transition clusters 6 and 8, group the MGnD microglia 7, 9, 11
Microglia_subset@meta.data$grouped_clusters <- case_when(Microglia_subset@meta.data$Clusters_microglia == "1" ~ "1",
                                                              Microglia_subset@meta.data$Clusters_microglia == "2" ~ "1",
                                                              Microglia_subset@meta.data$Clusters_microglia == "3" ~ "3",
                                                              Microglia_subset@meta.data$Clusters_microglia == "4" ~ "4",
                                                              Microglia_subset@meta.data$Clusters_microglia == "5" ~ "5",
                                                              Microglia_subset@meta.data$Clusters_microglia == "6" ~ "6",
                                                              Microglia_subset@meta.data$Clusters_microglia == "7" ~ "7",
                                                              Microglia_subset@meta.data$Clusters_microglia == "8" ~ "6",
                                                              Microglia_subset@meta.data$Clusters_microglia == "9" ~ "7",
                                                              Microglia_subset@meta.data$Clusters_microglia == "10" ~ "10",
                                                              Microglia_subset@meta.data$Clusters_microglia == "11" ~ "7", 
                                                              Microglia_subset@meta.data$Clusters_microglia == "12" ~ "12")

################## FEMALE AD dataset ####################################

Idents(Microglia_subset) <- "sex"
Microglia_female <- subset(x = Microglia_subset, subset = sex == "f")
table(Microglia_female@meta.data$group)
Idents(Microglia_subset) <- "diagnosis"
Microglia_female_AD <- subset(x = Microglia_female, subset = diagnosis == "AD")
table(Microglia_female_AD@meta.data$group)

######################## Grouped clusters #################################

# 1. Check effect of categorical variables on cluster proportions

######### Batch ############      
propeller_result <- propeller(clusters=Microglia_female_AD@meta.data$grouped_clusters, sample=Microglia_female_AD@meta.data$sample, group=Microglia_female_AD@meta.data$snRN.Aseq_batch)
fwrite(propeller_result , paste0(Output.dir.subset, "Propeller_results/Batch_female_AD.csv"), row.names = TRUE)
######### APOE ##############     
propeller_result <- propeller(clusters=Microglia_female_AD@meta.data$grouped_clusters, sample=Microglia_female_AD@meta.data$sample, group=Microglia_female_AD@meta.data$apoe)
fwrite(propeller_result , paste0(Output.dir.subset, "Propeller_results/APOE_genotype_female_AD.csv"), row.names = TRUE)

                                                         
# 2. Obtain clusters proportions 

# Obtain a dataframe with cluster proportions, counts, and transformed proportions per group

props <- getTransformedProps(clusters=Microglia_subset@meta.data$grouped_clusters, sample=Microglia_subset@meta.data$group, transform="logit")
fwrite(props$Proportions , paste0(Output.dir.subset, "Propeller_results/Cluster_proportions_whole_dataset_by_group_AD_female.csv"), row.names = TRUE)

# Plot proportions bar plot per group
pdf(file = paste0(Output.dir.subset, "Propeller_results/Cell_type_proportions_bar_plot_by_group.pdf"))
barplot(props$Proportions, col = c("orange","purple", "darkgreen", "red", "lightblue", "pink", "yellow", "grey", "lightgreen", "brown", "darkblue", "white", "gold"), 
        legend.text = TRUE,ylab = "Proportions", xlim=c(0,100))
dev.off()

# Obtain a dataframe with cluster proportions, counts, and transformed proportions per sample 

props <- getTransformedProps(clusters=Microglia_subset@meta.data$grouped_clusters, sample=Microglia_subset@meta.data$sample, transform="logit")
fwrite(props$Proportions , paste0(Output.dir.subset, "Propeller_results/Cluster_proportions_whole_dataset_AD_female.csv"), row.names = TRUE)

# Plot proportions bar plot per sample
pdf(file = paste0(Output.dir.subset, "Propeller_results/Cell_type_proportions_bar_plot_AD_female.pdf"))
barplot(props$Proportions, col = c("orange","purple", "darkgreen", "red", "lightblue", "pink", "yellow", "grey", "lightgreen", "brown", "darkblue", "white", "gold"), 
        legend.text = TRUE,ylab = "Proportions", xlim=c(0,50))
dev.off()





