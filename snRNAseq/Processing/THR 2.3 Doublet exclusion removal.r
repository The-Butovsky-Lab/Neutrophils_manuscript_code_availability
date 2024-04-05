# Donors have been genotyped for APOE, and have tau and amyloid pathology scores from the same tissue blocks used for snRNAseq. CellRanger count martices
# are stored in the input folder, with a seperate folder for each donor. The metadata .csv file is in the input folder. The output folder is created.

# 2.3 Run Scrublet to remove doublet nuclei from the dataset
# Author: Thomas Rust (THR) adapted from Mirjam Koster (MKO) and Astrid Alsema 
# Date 22/02/2024 (updated)

#Load packages
library(Seurat)

# Set output directory and read in the object that has the DoubletFinder doublets indicated but not removed 
Output.dir <- ("Output/")
total.list <- readRDS(file = paste0(Output.dir, "2. listed_samples doublets indicated.rds"))

# Load in Scrubblet's predictions for each sample
doublet_files <- list.files(path = paste0(Output.dir,"R_QC/Scrublet"), pattern = "predicted", full.names = T)

predicted_doublets <- list()
for(i in 1:length(doublet_files)){
  predicted_doublets[i] <- read.csv(file = doublet_files[i], header = F)
  # sanity check for conguence between seurat object and generated Scrublet scores
  if(ncol(total.list[[i]]) == length(predicted_doublets[[i]])){
    print(paste0("Both the Seurat object and the Scrubblet data have #Nuclei of ", ncol(total.list[[i]])))
  } else {
    print("There is an issue with the length of the predicted doublets")
  }
  
}

names(predicted_doublets) <- names(total.list)

# Add Scrubblet's prediction to the metadata of the Seurat object
pdf(file = paste0(Output.dir, "2.2 Doublet exclusion with Scrubblet.pdf"), width = 16, height = 9)
for(i in 1:length(total.list)){
  total.list[[i]] <- AddMetaData(object = total.list[[i]], metadata = predicted_doublets[[i]],
                                 col.name = "Scrubblet_predicted")

	
	print(DimPlot(total.list[[i]]) + DimPlot(total.list[[i]], group.by = "Scrubblet_predicted"))
    
}
dev.off()

# Save the object with the doublets indicated for both DoubletFinder and Scrublet scores, but no doublets removed
rm(list=setdiff(ls(), c("total.list", "Output.dir")))
saveRDS(object = total.list, file = paste0(Output.dir, "2.2 listed_samples doublets indicated.rds"))

## REMOVE THE DOUBLETS
for(i in 1:length(total.list)){
  total.list[[i]] <- subset(total.list[[i]], Scrubblet_predicted == 0)
}

# Save the seurat object with the doublets removed by Scrublet
saveRDS(total.list, file = paste0(Output.dir, "2.2 listed_samples doublets Scrublet removed.rds"))

# Proceed to THR 3.0 Integration of samples