# Donors have been genotyped for APOE, and have tau and amyloid pathology scores from the same tissue blocks used for snRNAseq. CellRanger count martices
# are stored in the input folder, with a seperate folder for each donor. The metadata3334 .csv file is in the input folder. The output folder is created. 
# 0. Load data and metadata - some metadata were not added at this stage but added later after object integration and annotation
# Author: Thomas Rust(THR) adpated from Mirjam Koster (MKO)
# Date: 14/01/24 (update)

# install old Seurat and SeuratObject versions to ensure version copatibility with pipeline
remotes::install_version("Seurat", version = "4.3.0")
remove.packages("SeuratObject") 
remotes::install_version("SeuratObject", version = "4.1.3")

# load packages 
library(Seurat)
library(reticulate)
library(dplyr)
library(patchwork)
library(ggplot2)
sessionInfo()

#Input directory contains folders named by sample name, containing the barcodes, features and matrix.mtx.gz files
Input.dir <- "Input"
Output.dir <- "Output/"

#############################################################################
### 0. Load data and metadata
#############################################################################

## LOADING THE CELLRANGER OUTPUT INTO SEURAT

# Load data
sample.input = list.dirs(path = Input.dir, recursive = F)
listed_sampledata <- lapply(X = sample.input, Read10X)

# Name the elements in the list based on folder names
sample.input.names = list.dirs(path = Input.dir, recursive = F, full.names = F)
names(listed_sampledata) <- sample.input.names

# Create Seurat objects from each element of the list
total.list <- lapply(X = listed_sampledata, FUN = function(x) {
  x <- CreateSeuratObject(x, min.cells = 3, min.features = 200, project = "snAD")
})

# Add a column in the metadata to identify the sample
for(a in 1:length(total.list)){
  total.list[[a]]$sample = names(total.list)[a]
  total.list[[a]] <- SetIdent(object = total.list[[a]], value = "sample")
}

## ADD METADATA FROM FILE INTO SEURAT OBJECTS

# Load metadata from a csv file
metadata0 <- read.csv2(file = paste0(Input.dir, "/metadata3334.csv"))
metadata0$sample <- as.character(metadata0$sample)
metadata0$Batch <- as.character(metadata0$Batch)
metadata0$Batch_Sample <- as.character(metadata0$Batch_Sample)

# Add the metadata to each individual sample in the list
for(i in 1:length(total.list)){
  metadata_total <- total.list[[i]]@meta.data[,c("orig.ident", "sample")]
  metadata_total$cell_id <- row.names(metadata_total)

  metadata <- right_join(x = metadata0, y = metadata_total, by = "sample")
  row.names(metadata) <- metadata$cell_id
  
  total.list[[i]] <- AddMetaData(object = total.list[[i]], metadata = metadata)
  
  write.table(x = total.list[[i]]@meta.data, sep = ";",
            file = paste0(Output.dir,"Metadata/0. Metadata sample ", names(total.list[i]),".csv"))
}

# Save
rm(list=setdiff(ls(), c("total.list", "Output.dir")))
saveRDS(object = total.list, file = paste0(Output.dir, "0. listed_samples unprocessed with metadata.rds"))