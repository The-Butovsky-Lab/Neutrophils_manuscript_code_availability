##
# Purpose: The purpose of the code is to perform FindAllMarkers on the Microglia object. 
# Date: February 20,2024
# Author: Madison Carpenter


# Setup the Enviornment
library(Seurat)
library(tidyverse)
library(magrittr)

# Assign Directories
exp.dir <- paste0("~/Dropbox (Partners HealthCare)/snRNAseq/Neta_APOE33_34_analysis/")
data.dir <- paste0("Output/6. Subcluster Microglia/")


# Load Data 
microglia <- readRDS(file = "~/Dropbox (Partners HealthCare)/snRNAseq/Neta_APOE33_34_analysis/Object/6.2 Microglia Round 3 resubclustered.rds")

## IDENTIFYING THE CLUSTERS - MAST
microglia.mast.cluster.markers <- FindAllMarkers(microglia, test.use = "MAST", latent.vars = c("sample"), assay = 'RNA', min.pct = 0.05)
microglia.mast.cluster.markers.DEGs <- microglia.mast.cluster.markers %>% filter(p_val <= 0.05)

writexl::write_xlsx(microglia.mast.cluster.markers,  "~/Dropbox (Partners HealthCare)/snRNAseq/Neta_APOE33_34_analysis/results/Microglia/FindAllMarkers/MG-APOE34-33_FindAllMarkers_MAST_nologthreshold.xlsx")
write.xlsx(microglia.mast.cluster.markers.DEGs,  "~/Dropbox (Partners HealthCare)/snRNAseq/Neta_APOE33_34_analysis/results/Microglia/FindAllMarkers/MG-APOE34-33_FindAllMarkers_MAST_nologthreshold_DEGs.xlsx")
