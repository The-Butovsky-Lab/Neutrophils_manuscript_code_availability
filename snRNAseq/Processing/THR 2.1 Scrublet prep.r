# Donors have been genotyped for APOE, and have tau and amyloid pathology scores from the same tissue blocks used for snRNAseq. CellRanger count martices
# are stored in the input folder, with a seperate folder for each donor. The metadata .csv file is in the input folder. The output folder is created.
# 2.1 Scrublet prep
# Author: Thomas Rust (THR) adapted from Mirjam Koster (MKO) and Astrid Alsema 
# Date 06-07-2023

library(Seurat)
library(MASS)
library(Matrix)
library(ggplot2)
library(stringr)
sessionInfo()

Output.dir <- "Output/"
total.list <- readRDS(file = paste0(Output.dir, "2. listed_samples doublets indicated.rds"))

### PREPARE DATA FOR SCRUBBLET (Python)

#Load the listed samples (preprocessed) with doublets from DoubletFinder indicated, but not removed
names(total.list)

dir.create(paste0(Output.dir,"R_QC"))
dir.create(paste0(Output.dir,"R_QC/Scrublet/"))
dir.create(paste0(Output.dir,"R_QC/Scrublet/counts"))

#Prepare matrices for loading into Scrubblet
for(i in 1:length(total.list)){
  dir.create(paste0(Output.dir, "R_QC/Scrublet/counts/", names(total.list)[i] ))
  writeMM(GetAssayData(total.list[[i]], slot = "counts"), paste(Output.dir, "R_QC/Scrublet/counts/", names(total.list)[i], "/matrix.mtx", sep = ""))
  barcodes <- Cells(total.list[[i]])
  features <- rownames(total.list[[i]]) 
  write.matrix(barcodes, paste(Output.dir, "R_QC/Scrublet/counts/", names(total.list)[i], "/barcodes.tsv", sep = ""))
  write.matrix(features, paste(Output.dir, "R_QC/Scrublet/counts/", names(total.list)[i], "/features.tsv", sep = ""))
}

# Proceed to THR 2.2 Scrublet score processing to run Scrublet