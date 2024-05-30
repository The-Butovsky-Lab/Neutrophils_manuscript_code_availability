########################################################################################
# Purpose: The purpose of the code is to perform FindAllMarkers on the Microglia object. 
# Date: February 20,2024
# Author: Madison Carpenter
########################################################################################

# Setup the Enviornment
library(Seurat)
library(tidyverse)
library(magrittr)

# Assign Directories
exp.dir <- paste0("Experiment Directory")
data.dir <- paste0("Output/6. Subcluster Microglia/")


# Load Data 
microglia <- readRDS(file = "Microglia Object RDS file")

## IDENTIFYING THE CLUSTERS - MAST
microglia.mast.cluster.markers <- FindAllMarkers(microglia, test.use = "MAST", latent.vars = c("sample"), assay = 'RNA', min.pct = 0.05)
microglia.mast.cluster.markers.DEGs <- microglia.mast.cluster.markers %>% filter(p_val <= 0.05)

writexl::write_xlsx(microglia.mast.cluster.markers,  "Microglia FindAllMarker results.xlsx")
write.xlsx(microglia.mast.cluster.markers.DEGs,  "Microglia FindAllMarkers Significant Results")
