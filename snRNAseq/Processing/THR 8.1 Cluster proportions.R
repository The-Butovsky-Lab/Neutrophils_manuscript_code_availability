# Author: Thomas Rust (THR)
# 8. Propeller for calculating cell type proportions
# Test metadata variables for effects on cell type proportions, then include the covariates in the final model
# Test these parameters in both the grouped and ungrouped methods for split male and female groups
# Date: 11/01/2024
################### Plotting known microglia markers ##################

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
Output.dir.subset = paste0("D:/Dropbox/snRNAseq/results/Microglia/CCA_Round3/Cell type proportions/APOE33_34/")
Microglia_subset <- readRDS(file = paste0("D:/Dropbox/snRNAseq/data/Microlgia/6.2 Microglia_Round3_metadata_apoe_dosage.rds"))

# Name the clusters
ID <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11") #### change number of clusters
names(ID) <- levels(Microglia_subset)
levels(Microglia_subset)
Microglia_subset$Clusters_microglia <- Idents(Microglia_subset)
Microglia_subset <- RenameIdents(Microglia_subset, ID)

## Assign cell categories for subclustering
# Group homeostatic clusters 0, 1 and 3 in to cluster 0 under a new column in the metadata - seurat_clusters_names
Microglia_subset@meta.data$seurat_clusters_names <- case_when(Microglia_subset@meta.data$Clusters_microglia == "0" ~ "0",
                                                              Microglia_subset@meta.data$Clusters_microglia == "1" ~ "0",
                                                              Microglia_subset@meta.data$Clusters_microglia == "2" ~ "2",
                                                              Microglia_subset@meta.data$Clusters_microglia == "3" ~ "0",
                                                              Microglia_subset@meta.data$Clusters_microglia == "4" ~ "4",
                                                              Microglia_subset@meta.data$Clusters_microglia == "5" ~ "5",
                                                              Microglia_subset@meta.data$Clusters_microglia == "6" ~ "6",
                                                              Microglia_subset@meta.data$Clusters_microglia == "7" ~ "7",
                                                              Microglia_subset@meta.data$Clusters_microglia == "8" ~ "8",
                                                              Microglia_subset@meta.data$Clusters_microglia == "9" ~ "9",
                                                              Microglia_subset@meta.data$Clusters_microglia == "10" ~ "10",
                                                              Microglia_subset@meta.data$Clusters_microglia == "11" ~ "11")

############# Subset data fro APOE33 and 34 donors ###########

Microglia_subset_APOE33_34 <- subset(Microglia_subset, apoe != "APOE44")

######################## Individual clusters - homeostatic clusters NOT grouped #################################

# 1. Check effect of categorical variables on cluster proportions

######### diagnosis ############   #### p < 0.05 (FDR) in clusters 5, 7 and 11 (almost cluster 4)
propeller_result <- propeller(clusters=Microglia_subset_APOE33_34@meta.data$Clusters_microglia, sample=Microglia_subset_APOE33_34@meta.data$sample, group=Microglia_subset_APOE33_34@meta.data$diagnosis)
fwrite(propeller_result , paste0(Output.dir.subset, "Individual clusters/Diagnosis_HM_not_grouped.csv"), row.names = TRUE)
######### Batch ############      ### p < 0.05 (FDR) in cluster 5 
propeller_result <- propeller(clusters=Microglia_subset_APOE33_34@meta.data$Clusters_microglia, sample=Microglia_subset_APOE33_34@meta.data$sample, group=Microglia_subset_APOE33_34@meta.data$snRNAseq_batch)
fwrite(propeller_result , paste0(Output.dir.subset, "Individual clusters/Batch_HM_not_grouped.csv"), row.names = TRUE)
######### apoe ##############     ### p > 0.05 (FDR) for all clusters, p < 0.05 in cluster 6 (increased APOE4 carriers)
propeller_result <- propeller(clusters=Microglia_subset_APOE33_34@meta.data$Clusters_microglia, sample=Microglia_subset_APOE33_34@meta.data$sample, group=Microglia_subset_APOE33_34@meta.data$apoe)
fwrite(propeller_result , paste0(Output.dir.subset, "Individual clusters/apoe_HM_not_grouped.csv"), row.names = TRUE)
######### apoe_3v4##############     ### p > 0.05 (FDR) for all clusters, p < 0.05 in clusters 6 and 3 (increased APOE4 carriers)
propeller_result <- propeller(clusters=Microglia_subset_APOE33_34@meta.data$Clusters_microglia, sample=Microglia_subset_APOE33_34@meta.data$sample, group=Microglia_subset_APOE33_34@meta.data$apoe_3v4)
fwrite(propeller_result , paste0(Output.dir.subset, "Individual clusters/apoe_3v4_HM_not_grouped.csv"), row.names = TRUE)
######### sex##############     ### p > 0.05 (FDR) for all clusters, p < 0.05 in cluster 6 (increased in females)
propeller_result <- propeller(clusters=Microglia_subset_APOE33_34@meta.data$Clusters_microglia, sample=Microglia_subset_APOE33_34@meta.data$sample, group=Microglia_subset_APOE33_34@meta.data$sex)
fwrite(propeller_result , paste0(Output.dir.subset, "Individual clusters/sex_HM_not_grouped.csv"), row.names = TRUE)

# 2. Obtain clusters proportions 
# Obtain a dataframe with cluster proportions, counts, and transformed proportions
props <- getTransformedProps(clusters=Microglia_subset_APOE33_34@meta.data$Clusters_microglia, sample=Microglia_subset_APOE33_34@meta.data$sample, transform="logit")
pdf(file = paste0(Output.dir.subset, "Individual clusters/Cell_type_proportions_bar_plot_HM_not_grouped.pdf"))
barplot(props$Proportions, col = c("orange","purple", "dark green", "red", "light blue", "pink", "yellow", "grey", "light green", "brown", "dark blue", "white", "gold"), 
        legend.text = TRUE,ylab = "Proportions", xlim=c(0,50))
dev.off()

# Read in sample info for design matrix to test effect of continous variables and construct designs with covariates included
sampleinfo <- read.csv2("D:/Dropbox/snRNAseq/Metadata/metadata_clust_apoe33_34.csv")
celltypes <- levels(Microglia_subset_APOE33_34)
group.microglia <- paste(sampleinfo$diagnosis, sampleinfo$pmd_min, sampleinfo$snRNAseq_batch, sep=",")

# Create metadata vectors as factors or numeric variables
grp <- as.factor(sampleinfo$diagnosis)
sex <- as.factor(sampleinfo$sex)
apoe <- as.factor(sampleinfo$apoe)
apoe3v4 <- as.factor(sampleinfo$apoe_3v4)
batch <- as.factor(sampleinfo$snRNAseq_batch)
pmd <- as.numeric(sampleinfo$pmd_min)
age <- as.numeric(sampleinfo$age)
tau <- as.numeric(sampleinfo$Tau_score)
ab <- as.numeric(sampleinfo$AB_score) # amyloid-beta has NAs, will need to test it manually or adjust model.matrix formula
prop.logit <- props
sampleinfo

## 3. Test effect of continous variables - we already know there is a batch effect

########## pmd ################   ### p > 0.05 (FDR) in all clusters
des.pmd <- model.matrix(~pmd)
fit <- lmFit(prop.logit$TransformedProps,des.pmd)
fit <- eBayes(fit, robust=TRUE)
eBayes_result <- topTable(fit,sort="none",n=Inf,coef="pmd")
fwrite(eBayes_result , paste0(Output.dir.subset, "Individual clusters/pmd_HM_not_grouped.csv"), row.names = TRUE)
######### age ################    ### p > 0.05 (FDR) in all clusters
des.age <- model.matrix(~age)
fit <- lmFit(prop.logit$TransformedProps,des.age)
fit <- eBayes(fit, robust=TRUE)
eBayes_result <- topTable(fit,fit,sort="none",n=Inf,coef="age")
fwrite(eBayes_result , paste0(Output.dir.subset, "Individual clusters/age_HM_not_grouped.csv"), row.names = TRUE)
##### no significant effect of age or pmd at FDR
######### age ################    ### p > 0.05 (FDR) in all clusters
des.tau <- model.matrix(~tau)
fit <- lmFit(prop.logit$TransformedProps,des.tau)
fit <- eBayes(fit, robust=TRUE)
eBayes_result <- topTable(fit,fit,sort="none",n=Inf,coef="tau")
fwrite(eBayes_result , paste0(Output.dir.subset, "Individual clusters/Tau_HM_not_grouped.csv"), row.names = TRUE)
##### 

###### 4. Test final cell proportion model 

### The only metadata parameter that was significantly influenced microglia subtype proportions was diagnosis
#Final proportioanl design - HOM clusters not grouped
design.anova <- model.matrix(~0+grp)
HM_not_grouped <- propeller.anova(prop.logit,design = design.anova, coef=c(1,2), robust=TRUE, 
                                  trend = FALSE, sort=TRUE)
fwrite(HM_not_grouped, paste0(Output.dir.subset, "Individual clusters/Propeller_final_diagnosis_HM_not_grouped.csv"), row.names = TRUE)

############### APOE as variable of interest - control for diagnosis and pmd

#Final proportioanl design - HOM clusters not grouped
design.anova <- model.matrix(~0+apoe+grp)
HM_not_grouped <- propeller.anova(prop.logit,design = design.anova, coef=c(1,2), robust=TRUE, 
                                  trend = FALSE, sort=TRUE)
fwrite(HM_not_grouped, paste0(Output.dir.subset, "Individual clusters/Propeller_APOE_final_HM_not_grouped.csv"), row.names = TRUE)


############# Subset data fro APOE33 and 34 donors ###########

Microglia_subset_APOE33_34 <- subset(Microglia_subset, apoe != "APOE44")

######################## Homeostatic clusters grouped #################################

# 1. Check effect of categorical variables on cluster proportions

######### diagnosis ############   #### p < 0.05 (FDR) in clusters 0, 4, 5, 7 and 11 
propeller_result <- propeller(clusters=Microglia_subset_APOE33_34@meta.data$seurat_clusters_names, sample=Microglia_subset_APOE33_34@meta.data$sample, group=Microglia_subset_APOE33_34@meta.data$diagnosis)
fwrite(propeller_result , paste0(Output.dir.subset, "Homeostatic clusters grouped/Diagnosis_HM_grouped.csv"), row.names = TRUE)
######### Batch ############      ### p < 0.05 (FDR) in cluster 5 
propeller_result <- propeller(clusters=Microglia_subset_APOE33_34@meta.data$seurat_clusters_names, sample=Microglia_subset_APOE33_34@meta.data$sample, group=Microglia_subset_APOE33_34@meta.data$snRNAseq_batch)
fwrite(propeller_result , paste0(Output.dir.subset, "Homeostatic clusters grouped/Batch_HM_grouped.csv"), row.names = TRUE)
######### apoe ##############     ### p > 0.05 (FDR) for all clusters, p < 0.05 in cluster 6 (increased APOE4 carriers)
propeller_result <- propeller(clusters=Microglia_subset_APOE33_34@meta.data$seurat_clusters_names, sample=Microglia_subset_APOE33_34@meta.data$sample, group=Microglia_subset_APOE33_34@meta.data$apoe)
fwrite(propeller_result , paste0(Output.dir.subset, "Homeostatic clusters grouped/apoe_HM_grouped.csv"), row.names = TRUE)
######### apoe_3v4##############     ### p > 0.05 (FDR) for all clusters, p < 0.05 in cluster 6 (increased APOE4 carriers)
propeller_result <- propeller(clusters=Microglia_subset_APOE33_34@meta.data$seurat_clusters_names, sample=Microglia_subset_APOE33_34@meta.data$sample, group=Microglia_subset_APOE33_34@meta.data$apoe_3v4)
fwrite(propeller_result , paste0(Output.dir.subset, "Homeostatic clusters grouped/apoe_3v4_HM_grouped.csv"), row.names = TRUE)
######### sex##############     ### p > 0.05 (FDR) for all clusters, p < 0.05 in cluster 6 (increased in females)
propeller_result <- propeller(clusters=Microglia_subset_APOE33_34@meta.data$seurat_clusters_names, sample=Microglia_subset_APOE33_34@meta.data$sample, group=Microglia_subset_APOE33_34@meta.data$sex)
fwrite(propeller_result , paste0(Output.dir.subset, "Homeostatic clusters grouped/sex_HM_grouped.csv"), row.names = TRUE)

# 2. Obtain clusters proportions 
# Obtain a dataframe with cluster proportions, counts, and transformed proportions
props <- getTransformedProps(clusters=Microglia_subset_APOE33_34@meta.data$seurat_clusters_names, sample=Microglia_subset_APOE33_34@meta.data$sample, transform="logit")
pdf(file = paste0(Output.dir.subset, "Homeostatic clusters grouped/Cell_type_proportions_bar_plot_HM_grouped.pdf"))
barplot(props$Proportions, col = c("orange","purple", "dark green", "red", "light blue", "pink", "yellow", "grey", "light green", "brown", "dark blue", "white", "gold"), 
        legend.text = TRUE,ylab = "Proportions", xlim=c(0,50))
dev.off()

# Read in sample info for design matrix to test effect of continous variables and construct designs with covariates included
sampleinfo <- read.csv2("D:/Dropbox/snRNAseq/Metadata/metadata_clust_apoe33_34.csv")
celltypes <- levels(Microglia_subset_APOE33_34)
group.microglia <- paste(sampleinfo$diagnosis, sampleinfo$pmd_min, sampleinfo$snRNAseq_batch, sep=",")

# Create metadata vectors as factors or numeric variables
grp <- as.factor(sampleinfo$diagnosis)
sex <- as.factor(sampleinfo$sex)
apoe <- as.factor(sampleinfo$apoe)
apoe3v4 <- as.factor(sampleinfo$apoe_3v4)
batch <- as.factor(sampleinfo$snRNAseq_batch)
pmd <- as.numeric(sampleinfo$pmd_min)
age <- as.numeric(sampleinfo$age)
tau <- as.numeric(sampleinfo$Tau_score)
ab <- as.numeric(sampleinfo$AB_score) # amyloid-beta has NAs, will need to test it manually or adjust model.matrix formula
prop.logit <- props
sampleinfo

## 3. Test effect of continous variables - we already know there is a batch effect

########## pmd ################   ### p > 0.05 (FDR) in all clusters
des.pmd <- model.matrix(~pmd)
fit <- lmFit(prop.logit$TransformedProps,des.pmd)
fit <- eBayes(fit, robust=TRUE)
eBayes_result <- topTable(fit,sort="none",n=Inf,coef="pmd")
fwrite(eBayes_result , paste0(Output.dir.subset, "Homeostatic clusters grouped/pmd_HM_grouped.csv"), row.names = TRUE)
######### age ################    ### p > 0.05 (FDR) in all clusters
des.age <- model.matrix(~age)
fit <- lmFit(prop.logit$TransformedProps,des.age)
fit <- eBayes(fit, robust=TRUE)
eBayes_result <- topTable(fit,fit,sort="none",n=Inf,coef="age")
fwrite(eBayes_result , paste0(Output.dir.subset, "Homeostatic clusters grouped/age_HM_grouped.csv"), row.names = TRUE)
##### no significant effect of age or pmd at FDR
######### age ################    ### p > 0.05 (FDR) in all clusters
des.tau <- model.matrix(~tau)
fit <- lmFit(prop.logit$TransformedProps,des.tau)
fit <- eBayes(fit, robust=TRUE)
eBayes_result <- topTable(fit,fit,sort="none",n=Inf,coef="tau")
fwrite(eBayes_result , paste0(Output.dir.subset, "Homeostatic clusters grouped/Tau_HM_grouped.csv"), row.names = TRUE)
##### 

###### 4. Test final cell proportion model 

### The only metadata parameter that was significantly influenced microglia subtype proportions was diagnosis 

#Final proportioanl design - HOM clusters not grouped
design.anova <- model.matrix(~0+grp)
HM_not_grouped <- propeller.anova(prop.logit,design = design.anova, coef=c(1,2), robust=TRUE, 
                                  trend = FALSE, sort=TRUE)
fwrite(HM_not_grouped, paste0(Output.dir.subset, "Homeostatic clusters grouped/Propeller_final_diagnosis_HM_grouped.csv"), row.names = TRUE)

############### APOE as variable of interest - control for diagnosis and batch

#Final proportioanl design - HOM clusters not grouped
design.anova <- model.matrix(~0+apoe+grp)
HM_not_grouped <- propeller.anova(prop.logit,design = design.anova, coef=c(1,2), robust=TRUE, 
                                  trend = FALSE, sort=TRUE)
fwrite(HM_not_grouped, paste0(Output.dir.subset, "Homeostatic clusters grouped/Propeller_APOE_final_HM_grouped.csv"), row.names = TRUE)