# Author: Thomas Rust (THR)
# 7.1. Propeller for calculating cell type proportions
# Test metadata variables for effects on cell type proportions, then include the covariates in the final model
# Test these parameters in individual clusters
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
Output.dir.subset = paste0("D:/Dropbox/snRNAseq/Neta_APOE33_34_analysis/results/Microglia/Cluster_proportions/")
Microglia_subset <- readRDS(file = paste0("D:/Dropbox/snRNAseq/Neta_APOE33_34_analysis/Object/6.2 MicrogliaFinalresubclustered_updated-ClusterNames.rds"))

# Name the clusters
ID <- c("1", "3", "8", "10", "4", "9", "2", "7", "5", "6", "11", "12") #### change number of clusters
names(ID) <- levels(Microglia_subset)
levels(Microglia_subset)
Microglia_subset$reordered_mappings <- Idents(Microglia_subset)
Microglia_subset <- RenameIdents(Microglia_subset, ID)
levels(Microglia_subset)


################## Whole dataset ####################################

######################## Individual clusters #################################

# 1. Check effect of categorical variables on cluster proportions

######### Diagnosis ############   
propeller_result <- propeller(clusters=Microglia_subset@meta.data$reordered_mappings, sample=Microglia_subset@meta.data$sample, group=Microglia_subset@meta.data$diagnosis)
fwrite(propeller_result , paste0(Output.dir.subset, "Individual clusters/Propeller_results/Diagnosis.csv"), row.names = TRUE)
######### Batch ############      
propeller_result <- propeller(clusters=Microglia_subset@meta.data$reordered_mappings, sample=Microglia_subset@meta.data$sample, group=Microglia_subset@meta.data$snRN.Aseq_batch)
fwrite(propeller_result , paste0(Output.dir.subset, "Individual clusters/Propeller_results/Batch.csv"), row.names = TRUE)
######### APOE ##############     
propeller_result <- propeller(clusters=Microglia_subset@meta.data$reordered_mappings, sample=Microglia_subset@meta.data$sample, group=Microglia_subset@meta.data$apoe)
fwrite(propeller_result , paste0(Output.dir.subset, "Individual clusters/Propeller_results/APOE_genotype.csv"), row.names = TRUE)
######### Sex ##############      
propeller_result <- propeller(clusters=Microglia_subset@meta.data$reordered_mappings, sample=Microglia_subset@meta.data$sample, group=Microglia_subset@meta.data$sex)
fwrite(propeller_result , paste0(Output.dir.subset, "Individual clusters/Propeller_results/Sex.csv"), row.names = TRUE)
######### Group ##############      
propeller_result <- propeller(clusters=Microglia_subset@meta.data$reordered_mappings, sample=Microglia_subset@meta.data$sample, group=Microglia_subset@meta.data$group)
fwrite(propeller_result , paste0(Output.dir.subset, "Individual clusters/Propeller_results/Group.csv"), row.names = TRUE)

# 2. Obtain clusters proportions 

# Obtain a dataframe with cluster proportions, counts, and transformed proportions per group

props <- getTransformedProps(clusters=Microglia_subset@meta.data$reordered_mappings, sample=Microglia_subset@meta.data$group, transform="logit")
fwrite(props$Proportions , paste0(Output.dir.subset, "Individual clusters/Propeller_results/Cluster_proportions_whole_dataset_by_group.csv"), row.names = TRUE)

# Plot proportions bar plot per group
pdf(file = paste0(Output.dir.subset, "Individual clusters/Propeller_results/Cell_type_proportions_bar_plot_by_group.pdf"))
barplot(props$Proportions, col = c("orange","purple", "darkgreen", "red", "lightblue", "pink", "yellow", "grey", "lightgreen", "brown", "darkblue", "white", "gold"), 
        legend.text = TRUE,ylab = "Proportions", xlim=c(0,100))
dev.off()

# Obtain a dataframe with cluster proportions, counts, and transformed proportions per sample 

props <- getTransformedProps(clusters=Microglia_subset@meta.data$reordered_mappings, sample=Microglia_subset@meta.data$sample, transform="logit")
fwrite(props$Proportions , paste0(Output.dir.subset, "Individual clusters/Propeller_results/Cluster_proportions_whole_dataset.csv"), row.names = TRUE)

# Plot proportions bar plot per sample
pdf(file = paste0(Output.dir.subset, "Individual clusters/Propeller_results/Cell_type_proportions_bar_plot.pdf"))
barplot(props$Proportions, col = c("orange","purple", "darkgreen", "red", "lightblue", "pink", "yellow", "grey", "lightgreen", "brown", "darkblue", "white", "gold"), 
        legend.text = TRUE,ylab = "Proportions", xlim=c(0,50))
dev.off()

# Read in sample info for design matrix to test effect of continous variables and construct designs with covariates included
sampleinfo <- read.csv2("D:/Dropbox/snRNAseq/Neta_APOE33_34_analysis/metadata/metadata_clust_apoe33_34.csv")
celltypes <- levels(Microglia_subset)
group.microglia <- paste(sampleinfo$diagnosis, sampleinfo$pmd_min, sampleinfo$snRNAseq_batch, sep=",")

# Create metadata vectors as factors or numeric variables
grp <- as.factor(sampleinfo$diagnosis)
sex <- as.factor(sampleinfo$sex)
apoe <- as.factor(sampleinfo$apoe)
batch <- as.factor(sampleinfo$snRNAseq_batch)
pmd <- as.numeric(sampleinfo$pmd_min)
age <- as.numeric(sampleinfo$age)
group <- as.numeric(sampleinfo$group)
tau <- as.numeric(sampleinfo$Tau_score)
ab <- as.numeric(sampleinfo$AB_score) # amyloid-beta has NAs, will need to test it manually or adjust model.matrix formula
prop.logit <- props
sampleinfo

## 3. Test effect of continous variables 

########## pmd ################   
des.pmd <- model.matrix(~pmd)
fit <- lmFit(prop.logit$TransformedProps,des.pmd)
fit <- eBayes(fit, robust=TRUE)
eBayes_result <- topTable(fit,sort="none",n=Inf,coef="pmd")
fwrite(eBayes_result , paste0(Output.dir.subset, "Individual clusters/Propeller_results/Post-mortem delay.csv"), row.names = TRUE)
######### age ################    
des.age <- model.matrix(~age)
fit <- lmFit(prop.logit$TransformedProps,des.age)
fit <- eBayes(fit, robust=TRUE)
eBayes_result <- topTable(fit,fit,sort="none",n=Inf,coef="age")
fwrite(eBayes_result , paste0(Output.dir.subset, "Individual clusters/Propeller_results/Age.csv"), row.names = TRUE)
##### no significant effect of age or pmd at FDR
######### age ################    
des.tau <- model.matrix(~tau)
fit <- lmFit(prop.logit$TransformedProps,des.tau)
fit <- eBayes(fit, robust=TRUE)
eBayes_result <- topTable(fit,fit,sort="none",n=Inf,coef="tau")
fwrite(eBayes_result , paste0(Output.dir.subset, "Individual clusters/Propeller_results/Tau.csv"), row.names = TRUE)


###### 4. Test final cell proportion model 

### The only metadata parameter that was significantly influenced microglia subtype proportions was diagnosis
#Final proportioanl design - HOM clusters not grouped
design.anova <- model.matrix(~0+grp)
model <- propeller.anova(prop.logit,design = design.anova, coef=c(1,2), robust=TRUE, 
                                  trend = FALSE, sort=TRUE)
fwrite(model, paste0(Output.dir.subset, "Individual clusters/Propeller_results/Propeller_final_diagnosis_effect.csv"), row.names = TRUE)

############### APOE as variable of interest - control for diagnosis and pmd

#Final proportioanl design - HOM clusters not grouped
design.anova <- model.matrix(~0+apoe+grp)
model <- propeller.anova(prop.logit,design = design.anova, coef=c(1,2), robust=TRUE, 
                                  trend = FALSE, sort=TRUE)
fwrite(model, paste0(Output.dir.subset, "Individual clusters/Propeller_results/Propeller_final_APOE_effect.csv"), row.names = TRUE)
