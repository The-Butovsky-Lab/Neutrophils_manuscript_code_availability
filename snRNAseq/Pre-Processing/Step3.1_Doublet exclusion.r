# Donors have been genotyped for APOE, and have tau and amyloid pathology scores from the same tissue blocks used for snRNAseq. 
# CellRanger count martices are stored in the input folder, with a seperate folder for each donor. The metadata .csv file is in the input folder. 
# The output folder is created.

# 2. Doublet exclusion with DoubletFinder. Doublet exlusion will also be performed with Scrublet. The efficacy of both methods can be compared and one 
# chosen to remove doublets. Scrublet seems to be working best and was used for this dataset

# Author: Thomas Rust adapted from Mirjam Koster (MKO)
# Date: 14/01/2024 (updated)


# Load packages
library(Seurat)
library(reticulate)
library(dplyr)
library(patchwork)
library(ggplot2)
library(DoubletFinder)

# Set output directory and read in processed object list of individual samples
Output.dir <- "Output/"
total.list <- readRDS(file = paste0(Output.dir, "1. listed_samples processed.rds"))

#############################################################################
### 2. Doublet exclusion
#############################################################################

pdf(file = paste0(Output.dir, "2. Doublet exclusion.pdf"), width = 16, height = 9)
for(i in 1:length(total.list)){
  
  #pK Identification
  sweep.res.list <- paramSweep_v3(total.list[[i]], PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  1
  chosen_pK = as.numeric(as.character(bcmvn$pK[bcmvn$BCmetric == max(bcmvn$BCmetric)]))
  
  print(ggplot(data = bcmvn, mapping = aes(pK, BCmetric, group = 1)) +
    geom_point() +
    geom_line() +
    theme_minimal() +
    geom_point(mapping = aes(x = pK[BCmetric == max(BCmetric)], y = max(BCmetric)), colour = "#49d22a") +
    ggtitle(paste0("sample ",names(total.list[i]),", chosen pK: ", chosen_pK)))

  #Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(annotations = total.list[[i]]@meta.data$seurat_clusters)
  nExp_poi <- round(0.075*nrow(total.list[[i]]@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies
  total.list[[i]] <- doubletFinder_v3(total.list[[i]], PCs = 1:30, pN = 0.25, pK = chosen_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  total.list[[i]] <- doubletFinder_v3(total.list[[i]], PCs = 1:30, pN = 0.25, pK = chosen_pK, nExp = nExp_poi.adj, reuse.pANN = paste0("pANN_0.25_", chosen_pK, "_", nExp_poi), sct = FALSE)
  
  total.list[[i]]$Doublet = paste0(total.list[[i]]@meta.data[,ncol(total.list[[i]]@meta.data)], "_",
                                   total.list[[i]]@meta.data[,ncol(total.list[[i]]@meta.data)-1])
  
  print(DimPlot(total.list[[i]]) + DimPlot(total.list[[i]], group.by = "Doublet"))
  
}
dev.off()

#Save the list in different forms.

## INCLUDE DOUBLETS - doublets marked but not removed
#rm(list=setdiff(ls(), c("total.list", "Output.dir")))
saveRDS(total.list, file = paste0(Output.dir, "2. listed_samples doublets indicated.rds"))

## REMOVE THE DOUBLETS STRINGENT - nuclei marked as doublet_doublet and Singlet_Doublet are both removed from the dataset
for(i in 1:length(total.list)){
  total.list[[i]] <- subset(total.list[[i]], Doublet == "Singlet_Singlet")
}
saveRDS(total.list, file = paste0(Output.dir, "2. listed_samples doublets removed stringent.rds"))

## REMOVE THE DOUBLETS - only nuclei marked as Doublet_Doublet are removed, therefore less stringent. Nuclei marked as Singlet_Doublet are included in this object
for(i in 1:length(total.list)){
  total.list[[i]] <- subset(total.list[[i]], Doublet != "Doublet_Doublet")
}
saveRDS(total.list, file = paste0(Output.dir, "2. listed_samples doublets removed.rds"))


# Proceed to using Scrublet with THR 2.1 Scrublet prep