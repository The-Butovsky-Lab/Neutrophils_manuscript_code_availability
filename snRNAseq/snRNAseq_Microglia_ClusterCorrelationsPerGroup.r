# 8. PLotting microglia correlation between cluster proportions and amyloid and tau pathology for certain groups
# Date: 29/02/2024

# LOAD PACKAGES
library(Seurat)
library(ggplot2)
library(clustree)
library(data.table)
library(readxl)
library(rstatix)
library(dplyr)
library(ggpubr)
library(heatmaply)
library(gplots)
library(binovisualfields)
library(corrplot)
library(ggpubr)
library(magrittr)
library(rstatix)
library(RColorBrewer)
sessionInfo()

# set information for cell type to analyse, set output directory and subdirectory
cell_type = "Microglia"
Round <- "Round 3"
Output.dir <- paste0("Output/")
Output.dir.subset = paste0(Output.dir, "6. Subcluster ",cell_type,"/Final/")
Output.dir.subset.results = paste0(Output.dir, "6. Subcluster ",cell_type,"/Final/")

# load in final microglia object
Microglia_subset <- readRDS(file = paste0(Output.dir.subset, "xxxx"))
Microglia_subset@meta.data$apoe <- as.factor(Microglia_subset@meta.data$apoe)

# Name the clusters
levels(Microglia_subset)

###############################################################################################################
###                                              Cluster distributions                                     ###
###############################################################################################################

microglia <- Microglia_subset

# calculate relative counts
rel_counts <- table(microglia@meta.data$reordered_mappings, microglia@meta.data$sample)
# calculate the relative counts
for (i in 1:ncol(rel_counts)){
  rel_counts[,i] <- (prop.table(rel_counts[,i]))*100 }
# make a dataframe of the relative counts
dg <- as.data.frame(rel_counts)
colnames(dg) <- c("cluster", "sample", "percentage")
# reorder factor levels
# to change the order on the x-axis of the barplot
dg$sample <- factor(dg$sample, levels = c("APOE10", "APOE20", "APOE30", "APOE40", "APOE50", "APOE60", "APOE19", "APOE29", "APOE39", 
                                          "APOE49", "APOE9", "APOE18", "APOE28", "APOE38", "APOE8", "APOE17", 
                                          "APOE27", "APOE37", "APOE47", "APOE57", "APOE7", "APOE3", "APOE13", "APOE23", 
                                          "APOE6", "APOE16", "APOE26", "APOE36", "APOE46", "APOE5", "APOE15", "APOE25", "APOE35", 
                                          "APOE45", "APOE4", "APOE14", "APOE24", "APOE34", "APOE44"))
write.csv(dg, file = paste0(Output.dir.subset.results, "7. Microglia cluster distribution ungrouped.csv"))
dg$sample <- as.factor(dg$sample)

# Add metadata to cluster distibution percentage table
metadatanew <- read.csv2("metadata/metadata_clust_apoe33_34.csv")

metadatanew$sample <- as.factor(metadatanew$sample)
metadatanew$apoe <- as.factor(metadatanew$apoe)
metadatanew$AB_score <- as.numeric(metadatanew$AB_score)
metadatanew$pmd_min <- as.numeric(metadatanew$pmd_min)
metadatanew$Tau_score <- as.numeric(metadatanew$Tau_score)
metadatanew$group <- as.factor(metadatanew$group)
metadatanew$group2 <- as.factor(metadatanew$group2)

# Add the metadata to each individual sample in the list                 
cell_clusters <- inner_join(dg, metadatanew, by = "sample")
cell_clusters$percentage <- as.numeric(cell_clusters$percentage)

### Subset clusters for analysis
cell_clusters$percentage <- as.numeric(cell_clusters$percentage)

micro1  <- subset(x = cell_clusters, subset = cluster == 1)
micro2  <- subset(x = cell_clusters, subset = cluster == 2)
micro3  <- subset(x = cell_clusters, subset = cluster == 3)
micro4  <- subset(x = cell_clusters, subset = cluster == 4)
micro5  <- subset(x = cell_clusters, subset = cluster == 5)
micro6  <- subset(x = cell_clusters, subset = cluster == 6)
micro7  <- subset(x = cell_clusters, subset = cluster == 7)
micro8  <- subset(x = cell_clusters, subset = cluster == 8)
micro9  <- subset(x = cell_clusters, subset = cluster == 9)
micro10  <- subset(x = cell_clusters, subset = cluster == 10)
micro11  <- subset(x = cell_clusters, subset = cluster == 11)
micro12  <- subset(x = cell_clusters, subset = cluster == 12)

# Define your color palette
color_palette <- brewer.pal(n = 12, "Set3")

# Define the mapping between micro objects and colors
micro_color_mapping <- c("micro1" = "#8DD3C7",
                         "micro2" = "#B3DE69",
                         "micro3" = "#FFFFB3",
                         "micro4" = "#80B1D3",
                         "micro5" = "#D9D9D9",
                         "micro6" = "#BC80BD",
                         "micro7" = "#FCCDE5",
                         "micro8" = "#BEBADA",
                         "micro9" = "#FDB462",
                         "micro10" = "#FB8072",
                         "micro11" = "#CCEBC5",
                         "micro12" =  "#FFED6F")

lineWidth <- 1
pointSize <- 20

# Tau correlation
for (i in 1:12) {
  micro_data <- get(paste0("micro", i))  # Get the data for micro object
  micro_name <- paste0("micro", i)  # Get the name of the micro object
  
  color_pal <- micro_color_mapping[micro_name]  # Select color from the mapping
  
  gg_corrlation <- ggplot(micro_data, aes(x = Tau_score, y = percentage)) +
    geom_point(size = 5, color = color_pal) +
    geom_smooth(method = "lm", se = TRUE, color = color_pal, formula = y ~ x) +
    stat_cor(color = "black", label.x.npc = "left", label.y.npc = "top") +
    theme(text = element_text(size = pointSize, colour = "black"),
          rect = element_blank(),
          line = element_line(size = lineWidth, colour = "black"),
          plot.title = element_text(color = "black", size = 20),
          axis.title.x = element_text(size = 15 , colour = "black"),
          axis.title.y = element_text(size = 15 , colour = "black"),
          axis.text.x = element_text(size = 15 , colour = "black", vjust = 0.1),
          axis.text.y = element_text(size = 15 , colour = "black"),
          axis.ticks = element_line(size = lineWidth, colour = "black"),
          axis.line = element_line(size = lineWidth, colour = "black"),
          legend.position = "right",
          legend.title = element_text(size = pointSize , colour = "black"),
          legend.text = element_text(size = pointSize , colour = "black"),
          legend.key.height = unit(1, "cm"),
          legend.key.width = unit(0.5, "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.margin = unit(0.5, "lines")) +
    ggtitle(paste("Cluster", i)) + 
    xlab("Tau score (%)") +
    ylab("Proportion of nuclei in cluster (%)") 
    
    
  # Save plot to PDF
  pdf(file = paste0(Output.dir.subset.results, "Micro_", i, "_Tau_correlation.pdf"))
  print(gg_corrlation)
  dev.off()
}

# AB correlation
for (i in 1:12) {
  micro_data <- get(paste0("micro", i))  # Get the data for micro object
  micro_name <- paste0("micro", i)  # Get the name of the micro object
  
  color_pal <- micro_color_mapping[micro_name]  # Select color from the mapping
  
  gg_corrlation <- ggplot(micro_data, aes(x = AB_score, y = percentage)) +
    geom_point(size = 5, color = color_pal) +
    geom_smooth(method = "lm", se = TRUE, color = color_pal, formula = y ~ x) +
    stat_cor(color = "black", label.x.npc = "left", label.y.npc = "top") +
    theme(text = element_text(size = pointSize, colour = "black"),
          rect = element_blank(),
          line = element_line(size = lineWidth, colour = "black"),
          plot.title = element_text(color = "black", size = 20),
          axis.title.x = element_text(size = 15 , colour = "black"),
          axis.title.y = element_text(size = 15 , colour = "black"),
          axis.text.x = element_text(size = 15 , colour = "black", vjust = 0.1),
          axis.text.y = element_text(size = 15 , colour = "black"),
          axis.ticks = element_line(size = lineWidth, colour = "black"),
          axis.line = element_line(size = lineWidth, colour = "black"),
          legend.position = "right",
          legend.title = element_text(size = pointSize , colour = "black"),
          legend.text = element_text(size = pointSize , colour = "black"),
          legend.key.height = unit(1, "cm"),
          legend.key.width = unit(0.5, "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.margin = unit(0.5, "lines")) +
    ggtitle(paste("Cluster", i)) + 
    xlab("AB score (%)") +
    ylab("Proportion of nuclei in cluster (%)") 
    
    
  # Save plot to PDF
  pdf(file = paste0(Output.dir.subset.results, "Micro_", i, "_AB_correlation.pdf"))
  print(gg_corrlation)
  dev.off()
}

# Subset females only
micro_female1 <- subset(x = micro1, subset = sex == "f")
micro_female2 <- subset(x = micro2, subset = sex == "f")
micro_female3 <- subset(x = micro3, subset = sex == "f")
micro_female4 <- subset(x = micro4, subset = sex == "f")
micro_female5 <- subset(x = micro5, subset = sex == "f")
micro_female6 <- subset(x = micro6, subset = sex == "f")
micro_female7 <- subset(x = micro7, subset = sex == "f")
micro_female8 <- subset(x = micro8, subset = sex == "f")
micro_female9 <- subset(x = micro9, subset = sex == "f")
micro_female10 <- subset(x = micro10, subset = sex == "f")
micro_female11 <- subset(x = micro11, subset = sex == "f")
micro_female12 <- subset(x = micro12, subset = sex == "f")

# Subset males only
micro_male1 <- subset(x = micro1, subset = sex == "m")
micro_male2 <- subset(x = micro2, subset = sex == "m")
micro_male3 <- subset(x = micro3, subset = sex == "m")
micro_male4 <- subset(x = micro4, subset = sex == "m")
micro_male5 <- subset(x = micro5, subset = sex == "m")
micro_male6 <- subset(x = micro6, subset = sex == "m")
micro_male7 <- subset(x = micro7, subset = sex == "m")
micro_male8 <- subset(x = micro8, subset = sex == "m")
micro_male9 <- subset(x = micro9, subset = sex == "m")
micro_male10 <- subset(x = micro10, subset = sex == "m")
micro_male11 <- subset(x = micro11, subset = sex == "m")
micro_male12 <- subset(x = micro12, subset = sex == "m")

########### Females ###################

# Define your color palette
color_palette <- brewer.pal(n = 12, "Set3")

# Define the mapping between micro objects and colors
micro_color_mapping <- c("micro_female1" = "#8DD3C7",
                         "micro_female2" = "#B3DE69",
                         "micro_female3" = "#FFFFB3",
                         "micro_female4" = "#80B1D3",
                         "micro_female5" = "#D9D9D9",
                         "micro_female6" = "#BC80BD",
                         "micro_female7" = "#FCCDE5",
                         "micro_female8" = "#BEBADA",
                         "micro_female9" = "#FDB462",
                         "micro_female10" = "#FB8072",
                         "micro_female11" = "#CCEBC5",
                         "micro_female12" =  "#FFED6F")

lineWidth <- 1
pointSize <- 20

# Tau correlation
for (i in 1:12) {
  micro_data <- get(paste0("micro_female", i))  # Get the data for micro object
  micro_name <- paste0("micro_female", i)  # Get the name of the micro object
  
  color_pal <- micro_color_mapping[micro_name]  # Select color from the mapping
  
  gg_corrlation <- ggplot(micro_data, aes(x = Tau_score, y = percentage)) +
    geom_point(size = 5, color = color_pal) +
    geom_smooth(method = "lm", se = TRUE, color = color_pal, formula = y ~ x) +
    stat_cor(color = "black", label.x.npc = "left", label.y.npc = "top") +
    theme(text = element_text(size = pointSize, colour = "black"),
          rect = element_blank(),
          line = element_line(size = lineWidth, colour = "black"),
          plot.title = element_text(color = "black", size = 20),
          axis.title.x = element_text(size = 15 , colour = "black"),
          axis.title.y = element_text(size = 15 , colour = "black"),
          axis.text.x = element_text(size = 15 , colour = "black", vjust = 0.1),
          axis.text.y = element_text(size = 15 , colour = "black"),
          axis.ticks = element_line(size = lineWidth, colour = "black"),
          axis.line = element_line(size = lineWidth, colour = "black"),
          legend.position = "right",
          legend.title = element_text(size = pointSize , colour = "black"),
          legend.text = element_text(size = pointSize , colour = "black"),
          legend.key.height = unit(1, "cm"),
          legend.key.width = unit(0.5, "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.margin = unit(0.5, "lines")) +
    ggtitle(paste("Cluster", i)) + 
    xlab("Tau score (%)") +
    ylab("Proportion of nuclei in cluster (%)") 
    
    
  # Save plot to PDF
  pdf(file = paste0(Output.dir.subset.results, "Micro_", i, "Female_Tau_correlation.pdf"))
  print(gg_corrlation)
  dev.off()
}

# AB correlation
for (i in 1:12) {
  micro_data <- get(paste0("micro_female", i))  # Get the data for micro object
  micro_name <- paste0("micro_female", i)  # Get the name of the micro object
  
  color_pal <- micro_color_mapping[micro_name]  # Select color from the mapping
  
  gg_corrlation <- ggplot(micro_data, aes(x = AB_score, y = percentage)) +
    geom_point(size = 5, color = color_pal) +
    geom_smooth(method = "lm", se = TRUE, color = color_pal, formula = y ~ x) +
    stat_cor(color = "black", label.x.npc = "left", label.y.npc = "top") +
    theme(text = element_text(size = pointSize, colour = "black"),
          rect = element_blank(),
          line = element_line(size = lineWidth, colour = "black"),
          plot.title = element_text(color = "black", size = 20),
          axis.title.x = element_text(size = 15 , colour = "black"),
          axis.title.y = element_text(size = 15 , colour = "black"),
          axis.text.x = element_text(size = 15 , colour = "black", vjust = 0.1),
          axis.text.y = element_text(size = 15 , colour = "black"),
          axis.ticks = element_line(size = lineWidth, colour = "black"),
          axis.line = element_line(size = lineWidth, colour = "black"),
          legend.position = "right",
          legend.title = element_text(size = pointSize , colour = "black"),
          legend.text = element_text(size = pointSize , colour = "black"),
          legend.key.height = unit(1, "cm"),
          legend.key.width = unit(0.5, "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.margin = unit(0.5, "lines")) +
    ggtitle(paste("Cluster", i)) + 
    xlab("AB score (%)") +
    ylab("Proportion of nuclei in cluster (%)") 
    
    
  # Save plot to PDF
  pdf(file = paste0(Output.dir.subset.results, "Micro_", i, "Female_AB_correlation.pdf"))
  print(gg_corrlation)
  dev.off()
}

########### Males ###################

# Define your color palette
color_palette <- brewer.pal(n = 12, "Set3")

# Define the mapping between micro objects and colors
micro_color_mapping <- c("micro_male1" = "#8DD3C7",
                         "micro_male2" = "#B3DE69",
                         "micro_male3" = "#FFFFB3",
                         "micro_male4" = "#80B1D3",
                         "micro_male5" = "#D9D9D9",
                         "micro_male6" = "#BC80BD",
                         "micro_male7" = "#FCCDE5",
                         "micro_male8" = "#BEBADA",
                         "micro_male9" = "#FDB462",
                         "micro_male10" = "#FB8072",
                         "micro_male11" = "#CCEBC5",
                         "micro_male12" =  "#FFED6F")

lineWidth <- 1
pointSize <- 20

# Tau correlation
for (i in 1:12) {
  micro_data <- get(paste0("micro_male", i))  # Get the data for micro object
  micro_name <- paste0("micro_male", i)  # Get the name of the micro object
  
  color_pal <- micro_color_mapping[micro_name]  # Select color from the mapping
  
  gg_corrlation <- ggplot(micro_data, aes(x = Tau_score, y = percentage)) +
    geom_point(size = 5, color = color_pal) +
    geom_smooth(method = "lm", se = TRUE, color = color_pal, formula = y ~ x) +
    stat_cor(color = "black", label.x.npc = "left", label.y.npc = "top") +
    theme(text = element_text(size = pointSize, colour = "black"),
          rect = element_blank(),
          line = element_line(size = lineWidth, colour = "black"),
          plot.title = element_text(color = "black", size = 20),
          axis.title.x = element_text(size = 15 , colour = "black"),
          axis.title.y = element_text(size = 15 , colour = "black"),
          axis.text.x = element_text(size = 15 , colour = "black", vjust = 0.1),
          axis.text.y = element_text(size = 15 , colour = "black"),
          axis.ticks = element_line(size = lineWidth, colour = "black"),
          axis.line = element_line(size = lineWidth, colour = "black"),
          legend.position = "right",
          legend.title = element_text(size = pointSize , colour = "black"),
          legend.text = element_text(size = pointSize , colour = "black"),
          legend.key.height = unit(1, "cm"),
          legend.key.width = unit(0.5, "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.margin = unit(0.5, "lines")) +
    ggtitle(paste("Cluster", i)) + 
    xlab("Tau score (%)") +
    ylab("Proportion of nuclei in cluster (%)") 
    
    
  # Save plot to PDF
  pdf(file = paste0(Output.dir.subset.results, "Micro_", i, "Male_Tau_correlation.pdf"))
  print(gg_corrlation)
  dev.off()
}

# AB correlation
for (i in 1:12) {
  micro_data <- get(paste0("micro_male", i))  # Get the data for micro object
  micro_name <- paste0("micro_male", i)  # Get the name of the micro object
  
  color_pal <- micro_color_mapping[micro_name]  # Select color from the mapping
  
  gg_corrlation <- ggplot(micro_data, aes(x = AB_score, y = percentage)) +
    geom_point(size = 5, color = color_pal) +
    geom_smooth(method = "lm", se = TRUE, color = color_pal, formula = y ~ x) +
    stat_cor(color = "black", label.x.npc = "left", label.y.npc = "top") +
    theme(text = element_text(size = pointSize, colour = "black"),
          rect = element_blank(),
          line = element_line(size = lineWidth, colour = "black"),
          plot.title = element_text(color = "black", size = 20),
          axis.title.x = element_text(size = 15 , colour = "black"),
          axis.title.y = element_text(size = 15 , colour = "black"),
          axis.text.x = element_text(size = 15 , colour = "black", vjust = 0.1),
          axis.text.y = element_text(size = 15 , colour = "black"),
          axis.ticks = element_line(size = lineWidth, colour = "black"),
          axis.line = element_line(size = lineWidth, colour = "black"),
          legend.position = "right",
          legend.title = element_text(size = pointSize , colour = "black"),
          legend.text = element_text(size = pointSize , colour = "black"),
          legend.key.height = unit(1, "cm"),
          legend.key.width = unit(0.5, "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.margin = unit(0.5, "lines")) +
    ggtitle(paste("Cluster", i)) + 
    xlab("AB score (%)") +
    ylab("Proportion of nuclei in cluster (%)") 
    
    
  # Save plot to PDF
  pdf(file = paste0(Output.dir.subset.results, "Micro_", i, "Male_AB_correlation.pdf"))
  print(gg_corrlation)
  dev.off()
}

################## FEMALES ALL CTRL AND AD33 ONLY ###############

# Subset females only - CTRL and AD33
micro_female_ctrl_AD33_1 <- subset(x = micro_female1, subset = group2 != "AD34")
micro_female_ctrl_AD33_2 <- subset(x = micro_female2, subset = group2 != "AD34")
micro_female_ctrl_AD33_3 <- subset(x = micro_female3, subset = group2 != "AD34")
micro_female_ctrl_AD33_4 <- subset(x = micro_female4, subset = group2 != "AD34")
micro_female_ctrl_AD33_5 <- subset(x = micro_female5, subset = group2 != "AD34")
micro_female_ctrl_AD33_6 <- subset(x = micro_female6, subset = group2 != "AD34")
micro_female_ctrl_AD33_7 <- subset(x = micro_female7, subset = group2 != "AD34")
micro_female_ctrl_AD33_8 <- subset(x = micro_female8, subset = group2 != "AD34")
micro_female_ctrl_AD33_9 <- subset(x = micro_female9, subset = group2 != "AD34")
micro_female_ctrl_AD33_10 <- subset(x = micro_female10, subset = group2 != "AD34")
micro_female_ctrl_AD33_11 <- subset(x = micro_female11, subset = group2 != "AD34")
micro_female_ctrl_AD33_12 <- subset(x = micro_female12, subset = group2 != "AD34")

########### Females ###################

# Define your color palette
color_palette <- brewer.pal(n = 12, "Set3")

# Define the mapping between micro objects and colors
micro_color_mapping <- c("micro_female_ctrl_AD33_1" = "#8DD3C7",
                         "micro_female_ctrl_AD33_2" = "#B3DE69",
                         "micro_female_ctrl_AD33_3" = "#FFFFB3",
                         "micro_female_ctrl_AD33_4" = "#80B1D3",
                         "micro_female_ctrl_AD33_5" = "#D9D9D9",
                         "micro_female_ctrl_AD33_6" = "#BC80BD",
                         "micro_female_ctrl_AD33_7" = "#FCCDE5",
                         "micro_female_ctrl_AD33_8" = "#BEBADA",
                         "micro_female_ctrl_AD33_9" = "#FDB462",
                         "micro_female_ctrl_AD33_10" = "#FB8072",
                         "micro_female_ctrl_AD33_11" = "#CCEBC5",
                         "micro_female_ctrl_AD33_12" =  "#FFED6F")

lineWidth <- 1
pointSize <- 20

# Tau correlation
for (i in 1:12) {
  micro_data <- get(paste0("micro_female_ctrl_AD33_", i))  # Get the data for micro object
  micro_name <- paste0("micro_female_ctrl_AD33_", i)  # Get the name of the micro object
  
  color_pal <- micro_color_mapping[micro_name]  # Select color from the mapping
  
  gg_corrlation <- ggplot(micro_data, aes(x = Tau_score, y = percentage)) +
    geom_point(size = 5, color = color_pal) +
    geom_smooth(method = "lm", se = TRUE, color = color_pal, formula = y ~ x) +
    stat_cor(color = "black", label.x.npc = "left", label.y.npc = "top") +
    theme(text = element_text(size = pointSize, colour = "black"),
          rect = element_blank(),
          line = element_line(size = lineWidth, colour = "black"),
          plot.title = element_text(color = "black", size = 20),
          axis.title.x = element_text(size = 15 , colour = "black"),
          axis.title.y = element_text(size = 15 , colour = "black"),
          axis.text.x = element_text(size = 15 , colour = "black", vjust = 0.1),
          axis.text.y = element_text(size = 15 , colour = "black"),
          axis.ticks = element_line(size = lineWidth, colour = "black"),
          axis.line = element_line(size = lineWidth, colour = "black"),
          legend.position = "right",
          legend.title = element_text(size = pointSize , colour = "black"),
          legend.text = element_text(size = pointSize , colour = "black"),
          legend.key.height = unit(1, "cm"),
          legend.key.width = unit(0.5, "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.margin = unit(0.5, "lines")) +
    ggtitle(paste("Cluster", i)) + 
    xlab("Tau score (%)") +
    ylab("Proportion of nuclei in cluster (%)") 
    
    
  # Save plot to PDF
  pdf(file = paste0(Output.dir.subset.results, "Micro_", i, "Female_ctrl_AD33_Tau_correlation.pdf"))
  print(gg_corrlation)
  dev.off()
}

# AB correlation
for (i in 1:12) {
  micro_data <- get(paste0("micro_female_ctrl_AD33_", i))  # Get the data for micro object
  micro_name <- paste0("micro_female_ctrl_AD33_", i)  # Get the name of the micro object
  
  color_pal <- micro_color_mapping[micro_name]  # Select color from the mapping
  
  gg_corrlation <- ggplot(micro_data, aes(x = AB_score, y = percentage)) +
    geom_point(size = 5, color = color_pal) +
    geom_smooth(method = "lm", se = TRUE, color = color_pal, formula = y ~ x) +
    stat_cor(color = "black", label.x.npc = "left", label.y.npc = "top") +
    theme(text = element_text(size = pointSize, colour = "black"),
          rect = element_blank(),
          line = element_line(size = lineWidth, colour = "black"),
          plot.title = element_text(color = "black", size = 20),
          axis.title.x = element_text(size = 15 , colour = "black"),
          axis.title.y = element_text(size = 15 , colour = "black"),
          axis.text.x = element_text(size = 15 , colour = "black", vjust = 0.1),
          axis.text.y = element_text(size = 15 , colour = "black"),
          axis.ticks = element_line(size = lineWidth, colour = "black"),
          axis.line = element_line(size = lineWidth, colour = "black"),
          legend.position = "right",
          legend.title = element_text(size = pointSize , colour = "black"),
          legend.text = element_text(size = pointSize , colour = "black"),
          legend.key.height = unit(1, "cm"),
          legend.key.width = unit(0.5, "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.margin = unit(0.5, "lines")) +
    ggtitle(paste("Cluster", i)) + 
    xlab("AB score (%)") +
    ylab("Proportion of nuclei in cluster (%)") 
    
    
  # Save plot to PDF
  pdf(file = paste0(Output.dir.subset.results, "Micro_", i, "Female_ctrl_AD33_AB_correlation.pdf"))
  print(gg_corrlation)
  dev.off()
}

################## FEMALES ALL CTRL AND AD34 ONLY ###############

# Subset females only - CTRL and AD34
micro_female_ctrl_AD34_1 <- subset(x = micro_female1, subset = group2 != "AD33")
micro_female_ctrl_AD34_2 <- subset(x = micro_female2, subset = group2 != "AD33")
micro_female_ctrl_AD34_3 <- subset(x = micro_female3, subset = group2 != "AD33")
micro_female_ctrl_AD34_4 <- subset(x = micro_female4, subset = group2 != "AD33")
micro_female_ctrl_AD34_5 <- subset(x = micro_female5, subset = group2 != "AD33")
micro_female_ctrl_AD34_6 <- subset(x = micro_female6, subset = group2 != "AD33")
micro_female_ctrl_AD34_7 <- subset(x = micro_female7, subset = group2 != "AD33")
micro_female_ctrl_AD34_8 <- subset(x = micro_female8, subset = group2 != "AD33")
micro_female_ctrl_AD34_9 <- subset(x = micro_female9, subset = group2 != "AD33")
micro_female_ctrl_AD34_10 <- subset(x = micro_female10, subset = group2 != "AD33")
micro_female_ctrl_AD34_11 <- subset(x = micro_female11, subset = group2 != "AD33")
micro_female_ctrl_AD34_12 <- subset(x = micro_female12, subset = group2 != "AD33")

########### Females ###################

# Define your color palette
color_palette <- brewer.pal(n = 12, "Set3")

# Define the mapping between micro objects and colors
micro_color_mapping <- c("micro_female_ctrl_AD34_1" = "#8DD3C7",
                         "micro_female_ctrl_AD34_2" = "#B3DE69",
                         "micro_female_ctrl_AD34_3" = "#FFFFB3",
                         "micro_female_ctrl_AD34_4" = "#80B1D3",
                         "micro_female_ctrl_AD34_5" = "#D9D9D9",
                         "micro_female_ctrl_AD34_6" = "#BC80BD",
                         "micro_female_ctrl_AD34_7" = "#FCCDE5",
                         "micro_female_ctrl_AD34_8" = "#BEBADA",
                         "micro_female_ctrl_AD34_9" = "#FDB462",
                         "micro_female_ctrl_AD34_10" = "#FB8072",
                         "micro_female_ctrl_AD34_11" = "#CCEBC5",
                         "micro_female_ctrl_AD34_12" =  "#FFED6F")

lineWidth <- 1
pointSize <- 20

# Tau correlation
for (i in 1:12) {
  micro_data <- get(paste0("micro_female_ctrl_AD34_", i))  # Get the data for micro object
  micro_name <- paste0("micro_female_ctrl_AD34_", i)  # Get the name of the micro object
  
  color_pal <- micro_color_mapping[micro_name]  # Select color from the mapping
  
  gg_corrlation <- ggplot(micro_data, aes(x = Tau_score, y = percentage)) +
    geom_point(size = 5, color = color_pal) +
    geom_smooth(method = "lm", se = TRUE, color = color_pal, formula = y ~ x) +
    stat_cor(color = "black", label.x.npc = "left", label.y.npc = "top") +
    theme(text = element_text(size = pointSize, colour = "black"),
          rect = element_blank(),
          line = element_line(size = lineWidth, colour = "black"),
          plot.title = element_text(color = "black", size = 20),
          axis.title.x = element_text(size = 15 , colour = "black"),
          axis.title.y = element_text(size = 15 , colour = "black"),
          axis.text.x = element_text(size = 15 , colour = "black", vjust = 0.1),
          axis.text.y = element_text(size = 15 , colour = "black"),
          axis.ticks = element_line(size = lineWidth, colour = "black"),
          axis.line = element_line(size = lineWidth, colour = "black"),
          legend.position = "right",
          legend.title = element_text(size = pointSize , colour = "black"),
          legend.text = element_text(size = pointSize , colour = "black"),
          legend.key.height = unit(1, "cm"),
          legend.key.width = unit(0.5, "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.margin = unit(0.5, "lines")) +
    ggtitle(paste("Cluster", i)) + 
    xlab("Tau score (%)") +
    ylab("Proportion of nuclei in cluster (%)") 
    
    
  # Save plot to PDF
  pdf(file = paste0(Output.dir.subset.results, "Micro_", i, "Female_ctrl_AD34_Tau_correlation.pdf"))
  print(gg_corrlation)
  dev.off()
}

# AB correlation
for (i in 1:12) {
  micro_data <- get(paste0("micro_female_ctrl_AD34_", i))  # Get the data for micro object
  micro_name <- paste0("micro_female_ctrl_AD34_", i)  # Get the name of the micro object
  
  color_pal <- micro_color_mapping[micro_name]  # Select color from the mapping
  
  gg_corrlation <- ggplot(micro_data, aes(x = AB_score, y = percentage)) +
    geom_point(size = 5, color = color_pal) +
    geom_smooth(method = "lm", se = TRUE, color = color_pal, formula = y ~ x) +
    stat_cor(color = "black", label.x.npc = "left", label.y.npc = "top") +
    theme(text = element_text(size = pointSize, colour = "black"),
          rect = element_blank(),
          line = element_line(size = lineWidth, colour = "black"),
          plot.title = element_text(color = "black", size = 20),
          axis.title.x = element_text(size = 15 , colour = "black"),
          axis.title.y = element_text(size = 15 , colour = "black"),
          axis.text.x = element_text(size = 15 , colour = "black", vjust = 0.1),
          axis.text.y = element_text(size = 15 , colour = "black"),
          axis.ticks = element_line(size = lineWidth, colour = "black"),
          axis.line = element_line(size = lineWidth, colour = "black"),
          legend.position = "right",
          legend.title = element_text(size = pointSize , colour = "black"),
          legend.text = element_text(size = pointSize , colour = "black"),
          legend.key.height = unit(1, "cm"),
          legend.key.width = unit(0.5, "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.margin = unit(0.5, "lines")) +
    ggtitle(paste("Cluster", i)) + 
    xlab("AB score (%)") +
    ylab("Proportion of nuclei in cluster (%)") 
    
    
  # Save plot to PDF
  pdf(file = paste0(Output.dir.subset.results, "Micro_", i, "Female_ctrl_AD34_AB_correlation.pdf"))
  print(gg_corrlation)
  dev.off()
}

################## MALES ALL CTRL AND AD33 ONLY ###############

# Subset males only - CTRL and AD33
micro_male_ctrl_AD33_1 <- subset(x = micro_male1, subset = group2 != "AD34")
micro_male_ctrl_AD33_2 <- subset(x = micro_male2, subset = group2 != "AD34")
micro_male_ctrl_AD33_3 <- subset(x = micro_male3, subset = group2 != "AD34")
micro_male_ctrl_AD33_4 <- subset(x = micro_male4, subset = group2 != "AD34")
micro_male_ctrl_AD33_5 <- subset(x = micro_male5, subset = group2 != "AD34")
micro_male_ctrl_AD33_6 <- subset(x = micro_male6, subset = group2 != "AD34")
micro_male_ctrl_AD33_7 <- subset(x = micro_male7, subset = group2 != "AD34")
micro_male_ctrl_AD33_8 <- subset(x = micro_male8, subset = group2 != "AD34")
micro_male_ctrl_AD33_9 <- subset(x = micro_male9, subset = group2 != "AD34")
micro_male_ctrl_AD33_10 <- subset(x = micro_male10, subset = group2 != "AD34")
micro_male_ctrl_AD33_11 <- subset(x = micro_male11, subset = group2 != "AD34")
micro_male_ctrl_AD33_12 <- subset(x = micro_male12, subset = group2 != "AD34")

########### Males ###################

# Define your color palette
color_palette <- brewer.pal(n = 12, "Set3")

# Define the mapping between micro objects and colors
micro_color_mapping <- c("micro_male_ctrl_AD33_1" = "#8DD3C7",
                         "micro_male_ctrl_AD33_2" = "#B3DE69",
                         "micro_male_ctrl_AD33_3" = "#FFFFB3",
                         "micro_male_ctrl_AD33_4" = "#80B1D3",
                         "micro_male_ctrl_AD33_5" = "#D9D9D9",
                         "micro_male_ctrl_AD33_6" = "#BC80BD",
                         "micro_male_ctrl_AD33_7" = "#FCCDE5",
                         "micro_male_ctrl_AD33_8" = "#BEBADA",
                         "micro_male_ctrl_AD33_9" = "#FDB462",
                         "micro_male_ctrl_AD33_10" = "#FB8072",
                         "micro_male_ctrl_AD33_11" = "#CCEBC5",
                         "micro_male_ctrl_AD33_12" =  "#FFED6F")

lineWidth <- 1
pointSize <- 20

# Tau correlation
for (i in 1:12) {
  micro_data <- get(paste0("micro_male_ctrl_AD33_", i))  # Get the data for micro object
  micro_name <- paste0("micro_male_ctrl_AD33_", i)  # Get the name of the micro object
  
  color_pal <- micro_color_mapping[micro_name]  # Select color from the mapping
  
  gg_corrlation <- ggplot(micro_data, aes(x = Tau_score, y = percentage)) +
    geom_point(size = 5, color = color_pal) +
    geom_smooth(method = "lm", se = TRUE, color = color_pal, formula = y ~ x) +
    stat_cor(color = "black", label.x.npc = "left", label.y.npc = "top") +
    theme(text = element_text(size = pointSize, colour = "black"),
          rect = element_blank(),
          line = element_line(size = lineWidth, colour = "black"),
          plot.title = element_text(color = "black", size = 20),
          axis.title.x = element_text(size = 15 , colour = "black"),
          axis.title.y = element_text(size = 15 , colour = "black"),
          axis.text.x = element_text(size = 15 , colour = "black", vjust = 0.1),
          axis.text.y = element_text(size = 15 , colour = "black"),
          axis.ticks = element_line(size = lineWidth, colour = "black"),
          axis.line = element_line(size = lineWidth, colour = "black"),
          legend.position = "right",
          legend.title = element_text(size = pointSize , colour = "black"),
          legend.text = element_text(size = pointSize , colour = "black"),
          legend.key.height = unit(1, "cm"),
          legend.key.width = unit(0.5, "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.margin = unit(0.5, "lines")) +
    ggtitle(paste("Cluster", i)) + 
    xlab("Tau score (%)") +
    ylab("Proportion of nuclei in cluster (%)") 
    
    
  # Save plot to PDF
  pdf(file = paste0(Output.dir.subset.results, "Micro_", i, "Male_ctrl_AD33_Tau_correlation.pdf"))
  print(gg_corrlation)
  dev.off()
}

# AB correlation
for (i in 1:12) {
  micro_data <- get(paste0("micro_male_ctrl_AD33_", i))  # Get the data for micro object
  micro_name <- paste0("micro_male_ctrl_AD33_", i)  # Get the name of the micro object
  
  color_pal <- micro_color_mapping[micro_name]  # Select color from the mapping
  
  gg_corrlation <- ggplot(micro_data, aes(x = AB_score, y = percentage)) +
    geom_point(size = 5, color = color_pal) +
    geom_smooth(method = "lm", se = TRUE, color = color_pal, formula = y ~ x) +
    stat_cor(color = "black", label.x.npc = "left", label.y.npc = "top") +
    theme(text = element_text(size = pointSize, colour = "black"),
          rect = element_blank(),
          line = element_line(size = lineWidth, colour = "black"),
          plot.title = element_text(color = "black", size = 20),
          axis.title.x = element_text(size = 15 , colour = "black"),
          axis.title.y = element_text(size = 15 , colour = "black"),
          axis.text.x = element_text(size = 15 , colour = "black", vjust = 0.1),
          axis.text.y = element_text(size = 15 , colour = "black"),
          axis.ticks = element_line(size = lineWidth, colour = "black"),
          axis.line = element_line(size = lineWidth, colour = "black"),
          legend.position = "right",
          legend.title = element_text(size = pointSize , colour = "black"),
          legend.text = element_text(size = pointSize , colour = "black"),
          legend.key.height = unit(1, "cm"),
          legend.key.width = unit(0.5, "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.margin = unit(0.5, "lines")) +
    ggtitle(paste("Cluster", i)) + 
    xlab("AB score (%)") +
    ylab("Proportion of nuclei in cluster (%)") 
    
    
  # Save plot to PDF
  pdf(file = paste0(Output.dir.subset.results, "Micro_", i, "Male_ctrl_AD33_AB_correlation.pdf"))
  print(gg_corrlation)
  dev.off()
}

################## MALES ALL CTRL AND AD34 ONLY ###############

# Subset males only - CTRL and AD34
micro_male_ctrl_AD34_1 <- subset(x = micro_male1, subset = group2 != "AD33")
micro_male_ctrl_AD34_2 <- subset(x = micro_male2, subset = group2 != "AD33")
micro_male_ctrl_AD34_3 <- subset(x = micro_male3, subset = group2 != "AD33")
micro_male_ctrl_AD34_4 <- subset(x = micro_male4, subset = group2 != "AD33")
micro_male_ctrl_AD34_5 <- subset(x = micro_male5, subset = group2 != "AD33")
micro_male_ctrl_AD34_6 <- subset(x = micro_male6, subset = group2 != "AD33")
micro_male_ctrl_AD34_7 <- subset(x = micro_male7, subset = group2 != "AD33")
micro_male_ctrl_AD34_8 <- subset(x = micro_male8, subset = group2 != "AD33")
micro_male_ctrl_AD34_9 <- subset(x = micro_male9, subset = group2 != "AD33")
micro_male_ctrl_AD34_10 <- subset(x = micro_male10, subset = group2 != "AD33")
micro_male_ctrl_AD34_11 <- subset(x = micro_male11, subset = group2 != "AD33")
micro_male_ctrl_AD34_12 <- subset(x = micro_male12, subset = group2 != "AD33")

########### Males ###################

# Define your color palette
color_palette <- brewer.pal(n = 12, "Set3")

# Define the mapping between micro objects and colors
micro_color_mapping <- c("micro_male_ctrl_AD34_1" = "#8DD3C7",
                         "micro_male_ctrl_AD34_2" = "#B3DE69",
                         "micro_male_ctrl_AD34_3" = "#FFFFB3",
                         "micro_male_ctrl_AD34_4" = "#80B1D3",
                         "micro_male_ctrl_AD34_5" = "#D9D9D9",
                         "micro_male_ctrl_AD34_6" = "#BC80BD",
                         "micro_male_ctrl_AD34_7" = "#FCCDE5",
                         "micro_male_ctrl_AD34_8" = "#BEBADA",
                         "micro_male_ctrl_AD34_9" = "#FDB462",
                         "micro_male_ctrl_AD34_10" = "#FB8072",
                         "micro_male_ctrl_AD34_11" = "#CCEBC5",
                         "micro_male_ctrl_AD34_12" =  "#FFED6F")

lineWidth <- 1
pointSize <- 20

# Tau correlation
for (i in 1:12) {
  micro_data <- get(paste0("micro_male_ctrl_AD34_", i))  # Get the data for micro object
  micro_name <- paste0("micro_male_ctrl_AD34_", i)  # Get the name of the micro object
  
  color_pal <- micro_color_mapping[micro_name]  # Select color from the mapping
  
  gg_corrlation <- ggplot(micro_data, aes(x = Tau_score, y = percentage)) +
    geom_point(size = 5, color = color_pal) +
    geom_smooth(method = "lm", se = TRUE, color = color_pal, formula = y ~ x) +
    stat_cor(color = "black", label.x.npc = "left", label.y.npc = "top") +
    theme(text = element_text(size = pointSize, colour = "black"),
          rect = element_blank(),
          line = element_line(size = lineWidth, colour = "black"),
          plot.title = element_text(color = "black", size = 20),
          axis.title.x = element_text(size = 15 , colour = "black"),
          axis.title.y = element_text(size = 15 , colour = "black"),
          axis.text.x = element_text(size = 15 , colour = "black", vjust = 0.1),
          axis.text.y = element_text(size = 15 , colour = "black"),
          axis.ticks = element_line(size = lineWidth, colour = "black"),
          axis.line = element_line(size = lineWidth, colour = "black"),
          legend.position = "right",
          legend.title = element_text(size = pointSize , colour = "black"),
          legend.text = element_text(size = pointSize , colour = "black"),
          legend.key.height = unit(1, "cm"),
          legend.key.width = unit(0.5, "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.margin = unit(0.5, "lines")) +
    ggtitle(paste("Cluster", i)) + 
    xlab("Tau score (%)") +
    ylab("Proportion of nuclei in cluster (%)") 
    
    
  # Save plot to PDF
  pdf(file = paste0(Output.dir.subset.results, "Micro_", i, "Male_ctrl_AD34_Tau_correlation.pdf"))
  print(gg_corrlation)
  dev.off()
}

# AB correlation
for (i in 1:12) {
  micro_data <- get(paste0("micro_male_ctrl_AD34_", i))  # Get the data for micro object
  micro_name <- paste0("micro_male_ctrl_AD34_", i)  # Get the name of the micro object
  
  color_pal <- micro_color_mapping[micro_name]  # Select color from the mapping
  
  gg_corrlation <- ggplot(micro_data, aes(x = AB_score, y = percentage)) +
    geom_point(size = 5, color = color_pal) +
    geom_smooth(method = "lm", se = TRUE, color = color_pal, formula = y ~ x) +
    stat_cor(color = "black", label.x.npc = "left", label.y.npc = "top") +
    theme(text = element_text(size = pointSize, colour = "black"),
          rect = element_blank(),
          line = element_line(size = lineWidth, colour = "black"),
          plot.title = element_text(color = "black", size = 20),
          axis.title.x = element_text(size = 15 , colour = "black"),
          axis.title.y = element_text(size = 15 , colour = "black"),
          axis.text.x = element_text(size = 15 , colour = "black", vjust = 0.1),
          axis.text.y = element_text(size = 15 , colour = "black"),
          axis.ticks = element_line(size = lineWidth, colour = "black"),
          axis.line = element_line(size = lineWidth, colour = "black"),
          legend.position = "right",
          legend.title = element_text(size = pointSize , colour = "black"),
          legend.text = element_text(size = pointSize , colour = "black"),
          legend.key.height = unit(1, "cm"),
          legend.key.width = unit(0.5, "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.margin = unit(0.5, "lines")) +
    ggtitle(paste("Cluster", i)) + 
    xlab("AB score (%)") +
    ylab("Proportion of nuclei in cluster (%)") 
    
    
  # Save plot to PDF
  pdf(file = paste0(Output.dir.subset.results, "Micro_", i, "Male_ctrl_AD34_AB_correlation.pdf"))
  print(gg_corrlation)
  dev.off()
}