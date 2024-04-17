# 7. PLotting microglia cluster proportions and correlation between cluster proportions and amyloid and tau pathology
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
Microglia_subset <- readRDS(file = paste0(Output.dir.subset, "6.2 MicrogliaFinalresubclustered_updated-ClusterNames.rds"))
Microglia_subset@meta.data$apoe <- as.factor(Microglia_subset@meta.data$apoe)

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
# Add the metadata to each individual sample in the list                 
cell_clusters <- inner_join(dg, metadatanew, by = "sample")
cell_clusters$percentage <- as.numeric(cell_clusters$percentage)

########################################################################
#			       FUNCTIONS 	               	       #
########################################################################
# barplot with the seurat colors per cluster
# have the percentages of each cluster in a dataframe
barp_perc <- function(dataframe){
  ggplot(data=dataframe, aes(x=sample, y= percentage, fill = cluster)) +
    geom_bar(stat="identity", position = position_fill(reverse = TRUE), width = 0.95)+
    scale_fill_manual(values = ucols)+
    xlab(NULL) +
    ylab("Cells") +
    scale_y_continuous(labels=scales::percent) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

BoxEmma <- function(cluster, title){
  Tmp <- tmp[tmp$cluster == cluster,]
  ggplot(Tmp, aes(x = diagnosis, y = percentage)) +
    geom_boxplot() +
    scale_colour_manual(values = c("blue", "red")) +
    geom_point() +
    theme_classic() +
    ggtitle(title)
}

### Subset clusters for analysis
cell_clusters$percentage <- as.numeric(cell_clusters$percentage)

# subset by cluster
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

#subset by cluster and sex - females
micro1_female <- subset(x = micro1, subset = sex == "f")
micro2_female <- subset(x = micro2, subset = sex == "f")
micro3_female <- subset(x = micro3, subset = sex == "f")
micro4_female <- subset(x = micro4, subset = sex == "f")
micro5_female <- subset(x = micro5, subset = sex == "f")
micro6_female <- subset(x = micro6, subset = sex == "f")
micro7_female <- subset(x = micro7, subset = sex == "f")
micro8_female <- subset(x = micro8, subset = sex == "f")
micro9_female <- subset(x = micro9, subset = sex == "f")
micro10_female <- subset(x = micro10, subset = sex == "f")
micro11_female <- subset(x = micro11, subset = sex == "f")
micro12_female <- subset(x = micro12, subset = sex == "f")

# subset by clusters - males
micro1_male <- subset(x = micro1, subset = sex == "m")
micro2_male <- subset(x = micro2, subset = sex == "m")
micro3_male <- subset(x = micro3, subset = sex == "m")
micro4_male <- subset(x = micro4, subset = sex == "m")
micro5_male <- subset(x = micro5, subset = sex == "m")
micro6_male <- subset(x = micro6, subset = sex == "m")
micro7_male <- subset(x = micro7, subset = sex == "m")
micro8_male <- subset(x = micro8, subset = sex == "m")
micro9_male <- subset(x = micro9, subset = sex == "m")
micro10_male <- subset(x = micro10, subset = sex == "m")
micro11_male <- subset(x = micro11, subset = sex == "m")
micro12_male <- subset(x = micro12, subset = sex == "m")

#################### Perform correlation tests between clusters and pathology ###########################

#################### Whole dataset ###############

# Create a vector to hold the micro data frames
micro_data <- list(micro1, micro2, micro3, micro4, micro5, micro6, micro7, micro8, micro9, micro10, micro11, micro12)

# Loop through micro datasets
for (i in 1:12) {
    # Perform correlation tests
    cor_test_result_tau <- cor.test(micro_data[[i]]$percentage, micro_data[[i]]$Tau_score, method = "pearson")
    cor_test_summary_tau <- data.frame(correlation_coefficient = cor_test_result_tau$estimate, 
                                       p_value = cor_test_result_tau$p.value,
                                       test_statistic = cor_test_result_tau$statistic)
    # Save results to CSV for Tau_score
    write.csv(cor_test_summary_tau , paste0(Output.dir.subset, "7. Micro", i, "_Tau.csv"), row.names = TRUE)

    cor_test_result_ab <- cor.test(micro_data[[i]]$percentage, micro_data[[i]]$AB_score, method = "pearson")
    cor_test_summary_ab <- data.frame(correlation_coefficient = cor_test_result_ab$estimate, 
                                      p_value = cor_test_result_ab$p.value,
                                      test_statistic = cor_test_result_ab$statistic)
    # Save results to CSV for AB_score
    write.csv(cor_test_summary_ab , paste0(Output.dir.subset, "7. Micro", i, "_AB.csv"), row.names = TRUE)
}

# Create a list to hold the lists of female and male objects
micro_female_list <- list(micro1_female, micro2_female, micro3_female, micro4_female, micro5_female, micro6_female, micro7_female, micro8_female, micro9_female, micro10_female, micro11_female, micro12_female)
micro_male_list <- list(micro1_male, micro2_male, micro3_male, micro4_male, micro5_male, micro6_male, micro7_male, micro8_male, micro9_male, micro10_male, micro11_male, micro12_male)

# Loop through female and male lists
for (i in 1:12) {
    # Perform correlation tests for female objects
    cor_test_result_tau_female <- cor.test(micro_female_list[[i]]$percentage, micro_female_list[[i]]$Tau_score, method = "pearson")
    cor_test_summary_tau_female <- data.frame(correlation_coefficient = cor_test_result_tau_female$estimate, 
                                       p_value = cor_test_result_tau_female$p.value,
                                       test_statistic = cor_test_result_tau_female$statistic)
    # Save results to CSV for Tau_score and female
    write.csv(cor_test_summary_tau_female, paste0(Output.dir.subset, "7. Micro", i, "_Female_Tau.csv"), row.names = TRUE)

    cor_test_result_ab_female <- cor.test(micro_female_list[[i]]$percentage, micro_female_list[[i]]$AB_score, method = "pearson")
    cor_test_summary_ab_female <- data.frame(correlation_coefficient = cor_test_result_ab_female$estimate, 
                                      p_value = cor_test_result_ab_female$p.value,
                                      test_statistic = cor_test_result_ab_female$statistic)
    # Save results to CSV for AB_score and female
    write.csv(cor_test_summary_ab_female, paste0(Output.dir.subset, "7. Micro", i, "_Female_AB.csv"), row.names = TRUE)

    # Perform correlation tests for male objects
    cor_test_result_tau_male <- cor.test(micro_male_list[[i]]$percentage, micro_male_list[[i]]$Tau_score, method = "pearson")
    cor_test_summary_tau_male <- data.frame(correlation_coefficient = cor_test_result_tau_male$estimate, 
                                       p_value = cor_test_result_tau_male$p.value,
                                       test_statistic = cor_test_result_tau_male$statistic)
    # Save results to CSV for Tau_score and male
    write.csv(cor_test_summary_tau_male, paste0(Output.dir.subset, "7. Micro", i, "_Male_Tau.csv"), row.names = TRUE)

    cor_test_result_ab_male <- cor.test(micro_male_list[[i]]$percentage, micro_male_list[[i]]$AB_score, method = "pearson")
    cor_test_summary_ab_male <- data.frame(correlation_coefficient = cor_test_result_ab_male$estimate, 
                                      p_value = cor_test_result_ab_male$p.value,
                                      test_statistic = cor_test_result_ab_male$statistic)
    # Save results to CSV for AB_score and male
    write.csv(cor_test_summary_ab_male, paste0(Output.dir.subset, "7. Micro", i, "_Male_AB.csv"), row.names = TRUE)
}

# Define your color palette
color_palette <- brewer.pal(n = 12, "Set3")

# Define the mapping between micro objects and colors
micro_color_mapping <- c("micro1 " = "#8DD3C7",
                         "micro2 " = "#B3DE69",
                         "micro3 " = "#FFFFB3",
                         "micro4 " = "#80B1D3",
                         "micro5 " = "#D9D9D9",
                         "micro6 " = "#BC80BD",
                         "micro7 " = "#FCCDE5",
                         "micro8 " = "#BEBADA",
                         "micro9 " = "#FDB462",
                         "micro10 " = "#FB8072",
                         "micro11 " = "#CCEBC5",
                         "micro12 " =  "#FFED6F")

# Iterate over each micro object

lineWidth <- 1
pointSize <- 20

# Iterate over each micro object
for (i in 1:12) {
  micro_data <- get(paste0("micro", i))  # Get the data for micro object
  micro_name <- paste0("micro", i)  # Get the name of the micro object
  
  color <- micro_color_mapping[micro_name]  # Select color from the mapping
  
  gg_corrlation <- ggplot(micro_data, aes(x = Tau_score, y = percentage)) +
    geom_point(size = 5, color = color) +
    geom_smooth(method = "lm", se = TRUE, color = paste0(color), formula = y ~ x) +
    stat_cor(color = paste0(color), label.x.npc = "left", label.y.npc = "top") +
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
    ggtitle("Cluster", i) + 
    xlab("Tau score (%)") +
    ylab("Proportion of nuclei in cluster (%)") 
    
  # Save plot to PDF
  pdf(file = paste0(Output.dir.subset.results, "Micro_", i, "Tau_correlation.pdf"))
  print(gg_corrlation)
  dev.off()
}

# Iterate over each micro object
for (i in 1:12) {
  micro_data <- get(paste0("micro", i))  # Get the data for micro object
  micro_name <- paste0("micro", i)  # Get the name of the micro object
  
  color_pal <- micro_color_mapping[micro_name]  # Select color from the mapping 
    
  gg_corrlation <- ggplot(micro_data, aes(x = AB_score, y = percentage, color = color_pal)) +
    geom_point(size = 5, color = color_pal) + scale_colour_identity() +
    geom_smooth(method = "lm", se = TRUE, color = paste0(color_pal), formula = y ~ x) + scale_colour_identity() +
    stat_cor(color = "black", label.x.npc = "left", label.y.npc = "top") +  scale_colour_identity() +
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
    ggtitle("Cluster", i) + 
    xlab("AB score (%)") +
    ylab("Proportion of nuclei in cluster (%)") 
    
  # Save plot to PDF
  pdf(file = paste0(Output.dir.subset.results, "Micro_", i, "_AB_correlation.pdf"))
  print(gg_corrlation)
  dev.off()
}

    ############## Micro1 ##############

# summary statistics
micro1 %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro1 %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro1 %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro1 %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc1 <- micro1 %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc1 <- pwc1 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 1_APOE.pdf"))
ggplot(micro1, aes(x=apoe, y=percentage, color = sex)) +
    geom_boxplot() + scale_colour_manual(values = c("#8DD3C7", "#8DD3C7")) +
    geom_point(color = "#8DD3C7") +
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
    ggtitle("Cluster 1") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc1) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

    ############## Micro1 female ##############

# summary statistics
micro1_female %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro1_female %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro1_female %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro1_female %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc1 <- micro1_female %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc1 <- pwc1 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 1_female_APOE.pdf"))
ggplot(micro1_female, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#8DD3C7") + scale_colour_identity() +
    geom_point(color = "#8DD3C7") +
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
    ggtitle("Cluster 1 - Females") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc1) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

    ############## Micro1 Male ##############

# summary statistics
micro1_male %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro1_male %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro1_male %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro1_male %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc1 <- micro1_male %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc1 <- pwc1 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 1_male_APOE.pdf"))
ggplot(micro1_male, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#8DD3C7") + scale_colour_identity() +
    geom_point(color = "#8DD3C7") +
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
    ggtitle("Cluster 1 - Males") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc1) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 2 ############

# summary statistics
micro2 %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro2 %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro2 %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro2 %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc2 <- micro2 %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc2 <- pwc2 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 2_APOE.pdf"))
ggplot(micro2, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#B3DE69") + scale_colour_identity() +
    geom_point(color = "#B3DE69") +
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
    ggtitle("Cluster 2") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc2) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

    ############## Micro2 female ##############

# summary statistics
micro2_female %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro2_female %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro2_female %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro2_female %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc2 <- micro2_female %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc2 <- pwc2 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 2_female_APOE.pdf"))
ggplot(micro2_female, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#B3DE69") + scale_colour_identity() +
    geom_point(color = "#B3DE69") +
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
    ggtitle("Cluster 1 - Females") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc2) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

    ############## Micro2 Male ##############

# summary statistics
micro2_male %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro2_male %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro2_male %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro2_male %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc2 <- micro2_male %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc2 <- pwc2 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 2_male_APOE.pdf"))
ggplot(micro2_male, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#B3DE69") + scale_colour_identity() +
    geom_point(color = "#B3DE69") +
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
    ggtitle("Cluster 1 - Males") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc2) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 3 ############

# summary statistics
micro3 %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro3 %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro3 %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro3 %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc3 <- micro3 %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc3 <- pwc3 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 3_APOE.pdf"))
ggplot(micro3, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#FFFFB3") + scale_colour_identity() +
    geom_point(color = "#FFFFB3") +
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
    ggtitle("Cluster 3") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc3) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 3 - female ############

# summary statistics
micro3_female %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro3_female %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro3_female %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro3_female %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc3 <- micro3_female %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc3 <- pwc3 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 3_female_APOE.pdf"))
ggplot(micro3_female, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#FFFFB3") + scale_colour_identity() +
    geom_point(color = "#FFFFB3") +
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
    ggtitle("Cluster 3 -  Females") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc3) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 3 - male ############

# summary statistics
micro3_male %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro3_male %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro3_male %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro3_male %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc3 <- micro3_male %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc3 <- pwc3 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 3_male_APOE.pdf"))
ggplot(micro3_female, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#FFFFB3") + scale_colour_identity() +
    geom_point(color = "#FFFFB3") +
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
    ggtitle("Cluster 3 -  Males") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc3) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 4 ############

# summary statistics
micro4 %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro4 %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro4 %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro4 %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc4 <- micro4 %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc4 <- pwc4 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 4_APOE.pdf"))
ggplot(micro4, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#80B1D3") + scale_colour_identity() +
    geom_point(color = "#80B1D3") +
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
    ggtitle("Cluster 4") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc4) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 4 female ############

# summary statistics
micro4_female %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro4_female %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro4_female %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro4_female %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc4 <- micro4_female %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc4 <- pwc4 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 4_female_APOE.pdf"))
ggplot(micro4_female, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#80B1D3") + scale_colour_identity() +
    geom_point(color = "#80B1D3") +
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
    ggtitle("Cluster 4 - Females") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc4) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 4 male ############

# summary statistics
micro4_male %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro4_male %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro4_male %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro4_male %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc4 <- micro4_male %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc4 <- pwc4 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 4_male_APOE.pdf"))
ggplot(micro4_male, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#80B1D3") + scale_colour_identity() +
    geom_point(color = "#80B1D3") +
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
    ggtitle("Cluster 4 - Males") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc4) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()


######### Micro 5 ############

# summary statistics
micro5 %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro5 %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro5 %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro5 %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc5 <- micro5 %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc5 <- pwc5 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 5_APOE.pdf"))
ggplot(micro5, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#D9D9D9") + scale_colour_identity() +
    geom_point(color = "#D9D9D9") +
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
    ggtitle("Cluster 5") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc5) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 5 - females ############

# summary statistics
micro5_female %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro5_female %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro5_female %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro5_female %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc5 <- micro5_female %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc5 <- pwc5 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 5_female_APOE.pdf"))
ggplot(micro5_female, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#D9D9D9") + scale_colour_identity() +
    geom_point(color = "#D9D9D9") +
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
    ggtitle("Cluster 5 - Females") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc5) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 5 - males ############

# summary statistics
micro5_male %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro5_male %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro5_male %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro5_male %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc5 <- micro5_male %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc5 <- pwc5 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 5_male_APOE.pdf"))
ggplot(micro5_male, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#D9D9D9") + scale_colour_identity() +
    geom_point(color = "#D9D9D9") +
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
    ggtitle("Cluster 5 - Males") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc5) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 6 ############

# summary statistics
micro6 %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro6 %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro6 %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro6 %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc6 <- micro6 %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc6 <- pwc6 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 6_APOE.pdf"))
ggplot(micro6, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#BC80BD") + scale_colour_identity() +
    geom_point(color = "#BC80BD") +
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
    ggtitle("Cluster 6") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc6) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 6 - females ############

# summary statistics
micro6_female %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro6_female %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro6_female %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro6_female %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc6 <- micro6_female %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc6 <- pwc6 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 6_female_APOE.pdf"))
ggplot(micro6_female, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#BC80BD") + scale_colour_identity() +
    geom_point(color = "#BC80BD") +
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
    ggtitle("Cluster 6 - Females") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc6) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 6 - males ############

# summary statistics
micro6_male %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro6_male %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro6_male %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro6_male %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc6 <- micro6_male %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc6 <- pwc6 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 6_male_APOE.pdf"))
ggplot(micro6_male, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#BC80BD") + scale_colour_identity() +
    geom_point(color = "#BC80BD") +
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
    ggtitle("Cluster 6 - Males") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc6) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 7 ############

# summary statistics
micro7 %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro7 %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro7 %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro7 %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc7 <- micro7 %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc7 <- pwc7 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 7_APOE.pdf"))
ggplot(micro7, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#FCCDE5") + scale_colour_identity() +
    geom_point(color = "#FCCDE5") +
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
    ggtitle("Cluster 7") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc7) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 7 - females ############

# summary statistics
micro7_female %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro7_female %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro7_female %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro7_female %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc7 <- micro7_female %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc7 <- pwc7 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 7_female_APOE.pdf"))
ggplot(micro7_female, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#FCCDE5") + scale_colour_identity() +
    geom_point(color = "#FCCDE5") +
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
    ggtitle("Cluster 7 - Females") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc7) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 7 - males ############

# summary statistics
micro7_male %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro7_male %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro7_male %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro7_male %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc7 <- micro7_male %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc7 <- pwc7 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 7_male_APOE.pdf"))
ggplot(micro7_male, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#FCCDE5") + scale_colour_identity() +
    geom_point(color = "#FCCDE5") +
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
    ggtitle("Cluster 7 - Males") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc7) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 8 ############

# summary statistics
micro8 %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro8 %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro8 %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro8 %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc8 <- micro8 %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc8 <- pwc8 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 8_APOE.pdf"))
ggplot(micro8, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#BEBADA") + scale_colour_identity() +
    geom_point(color = "#BEBADA") +
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
    ggtitle("Cluster 8") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc8) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 8 - females ############

# summary statistics
micro8_female %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro8_female %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro8_female %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro8_female %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc8 <- micro8_female %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc8 <- pwc8 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 8_female_APOE.pdf"))
ggplot(micro8_female, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#BEBADA") + scale_colour_identity() +
    geom_point(color = "#BEBADA") +
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
    ggtitle("Cluster 8 - Females") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc8) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 8 ############

# summary statistics
micro8_male %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro8_male %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro8_male %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro8_male %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc8 <- micro8_male %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc8 <- pwc8 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 8_male_APOE.pdf"))
ggplot(micro8_male, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#BEBADA") + scale_colour_identity() +
    geom_point(color = "#BEBADA") +
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
    ggtitle("Cluster 8 - Males") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc8) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 9 ############

# summary statistics
micro9 %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro9 %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro9 %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro9 %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc9 <- micro9 %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc9 <- pwc9 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 9_APOE.pdf"))
ggplot(micro9, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#FDB462") + scale_colour_identity() +
    geom_point(color = "#FDB462") +
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
    ggtitle("Cluster 9") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc9) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 9 - females ############

# summary statistics
micro9_female %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro9_female %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro9_female %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro9_female %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc9 <- micro9_female %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc9 <- pwc9 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 9_female_APOE.pdf"))
ggplot(micro9_female, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#FDB462") + scale_colour_identity() +
    geom_point(color = "#FDB462") +
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
    ggtitle("Cluster 9 - Females") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc9) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 9 - males############

# summary statistics
micro9_male %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro9_male %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro9_male %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro9_male %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc9 <- micro9_male %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc9 <- pwc9 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 9_male_APOE.pdf"))
ggplot(micro9_male, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#FDB462") + scale_colour_identity() +
    geom_point(color = "#FDB462") +
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
    ggtitle("Cluster 9 - Males") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc9) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 10 ############

# summary statistics
micro10 %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro10 %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro10 %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro10 %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc10 <- micro10 %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc10 <- pwc10 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 10_APOE.pdf"))
ggplot(micro10, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#FB8072") + scale_colour_identity() +
    geom_point(color = "#FB8072") +
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
    ggtitle("Cluster 10") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc10) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 10 - females ############

# summary statistics
micro10_female %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro10_female %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro10_female %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro10_female %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc10 <- micro10_female %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc10 <- pwc10 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 10_female_APOE.pdf"))
ggplot(micro10_female, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#FB8072") + scale_colour_identity() +
    geom_point(color = "#FB8072") +
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
    ggtitle("Cluster 10 - Females") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc10) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 10 - males ############

# summary statistics
micro10_male %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro10_male %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro10_male %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro10_male %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc10 <- micro10_male %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc10 <- pwc10 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 10_male_APOE.pdf"))
ggplot(micro10_male, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#FB8072") + scale_colour_identity() +
    geom_point(color = "#FB8072") +
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
    ggtitle("Cluster 10 - Males") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc10) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 11 ############

# summary statistics
micro11 %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro11 %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro11 %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro11 %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc11 <- micro11 %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc11 <- pwc11 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 11_APOE.pdf"))
ggplot(micro11, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#CCEBC5") + scale_colour_identity() +
    geom_point(color = "#CCEBC5") +
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
    ggtitle("Cluster 11") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc11) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 11 - females ############

# summary statistics
micro11_female %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro11_female %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro11_female %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro11_female %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc11 <- micro11_female %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc11 <- pwc11 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 11_female_APOE.pdf"))
ggplot(micro11_female, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#CCEBC5") + scale_colour_identity() +
    geom_point(color = "#CCEBC5") +
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
    ggtitle("Cluster 11 - Females") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc11) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 11 - males ############

# summary statistics
micro11_male %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro11_male %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro11_male %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro11_male %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc11 <- micro11_male %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc11 <- pwc11 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 11_male_APOE.pdf"))
ggplot(micro11_male, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#CCEBC5") + scale_colour_identity() +
    geom_point(color = "#CCEBC5") +
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
    ggtitle("Cluster 11 - Males") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc11) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 12 ############

# summary statistics
micro12 %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro12 %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro12 %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro12 %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc12 <- micro12 %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc12 <- pwc12 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 12_APOE.pdf"))
ggplot(micro12, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#FFED6F") + scale_colour_identity() +
    geom_point(color = "#FFED6F") +
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
    ggtitle("Cluster 12") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc12) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 12 - females ############

# summary statistics
micro12_female %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro12_female %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro12_female %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro12_female %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc12 <- micro12_female %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc12 <- pwc12 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 12_female_APOE.pdf"))
ggplot(micro12_female, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#FFED6F") + scale_colour_identity() +
    geom_point(color = "#FFED6F") +
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
    ggtitle("Cluster 12 - Females") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc12) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 12 - males ############

# summary statistics
micro12_male %>% group_by(apoe) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro12_male %>%
    group_by(apoe) %>%
    identify_outliers(percentage)
# normality assumption
  micro12_male %>%
    group_by(apoe) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro12_male %>% anova_test(percentage ~  apoe)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc12 <- micro12_male %>% tukey_hsd(percentage ~  apoe)

# plot results
  # Visualization: box plots with p-values
  pwc12 <- pwc12 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 12_male_APOE.pdf"))
ggplot(micro12_male, aes(x=apoe, y=percentage)) +
    geom_boxplot(color = "#FFED6F") + scale_colour_identity() +
    geom_point(color = "#FFED6F") +
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
    ggtitle("Cluster 12 - Males") + 
    xlab("APOE genotype") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc12) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

############### BOXPLOTS BY DIAGNOSIS #####################

    ############## Micro1 ##############

# summary statistics
micro1 %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro1 %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro1 %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro1 %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc1 <- micro1 %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc1 <- pwc1 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 1_diagnosis.pdf"))
ggplot(micro1, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#8DD3C7") + scale_colour_identity() +
    geom_point(color = "#8DD3C7") +
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
    ggtitle("Cluster 1") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc1) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

    ############## Micro1 - females ##############

# summary statistics
micro1_female %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro1_female %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro1_female %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro1_female %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc1 <- micro1_female %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc1 <- pwc1 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 1_female_diagnosis.pdf"))
ggplot(micro1_female, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#8DD3C7") + scale_colour_identity() +
    geom_point(color = "#8DD3C7") +
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
    ggtitle("Cluster 1 - Females") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc1) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

    ############## Micro1 - females ##############

# summary statistics
micro1_female %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro1_female %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro1_female %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro1_female %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc1 <- micro1_female %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc1 <- pwc1 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 1_female_diagnosis.pdf"))
ggplot(micro1_female, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#8DD3C7") + scale_colour_identity() +
    geom_point(color = "#8DD3C7") +
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
    ggtitle("Cluster 1 - Females") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc1) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

    ############## Micro1 - males ##############

# summary statistics
micro1_male %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro1_male %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro1_male %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro1_male %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc1 <- micro1_male %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc1 <- pwc1 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 1_male_diagnosis.pdf"))
ggplot(micro1_male, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#8DD3C7") + scale_colour_identity() +
    geom_point(color = "#8DD3C7") +
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
    ggtitle("Cluster 1 - Males") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc1) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 2 ############

# summary statistics
micro2 %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro2 %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro2 %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro2 %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc2 <- micro2 %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc2 <- pwc2 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 2_diagnosis.pdf"))
ggplot(micro2, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#B3DE69") + scale_colour_identity() +
    geom_point(color = "#B3DE69") +
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
    ggtitle("Cluster 2") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc2) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 2 - females ############

# summary statistics
micro2_female %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro2_female %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro2_female %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro2_female %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc2 <- micro2_female %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc2 <- pwc2 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 2_female_diagnosis.pdf"))
ggplot(micro2_female, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#B3DE69") + scale_colour_identity() +
    geom_point(color = "#B3DE69") +
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
    ggtitle("Cluster 2 - Females") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc2) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 2 - males ############

# summary statistics
micro2_male %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro2_male %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro2_male %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro2_male %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc2 <- micro2_male %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc2 <- pwc2 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 2_male_diagnosis.pdf"))
ggplot(micro2_male, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#B3DE69") + scale_colour_identity() +
    geom_point(color = "#B3DE69") +
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
    ggtitle("Cluster 2 - Males") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc2) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 3 ############

# summary statistics
micro3 %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro3 %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro3 %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro3 %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc3 <- micro3 %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc3 <- pwc3 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 3_diagnosis.pdf"))
ggplot(micro3, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#FFFFB3") + scale_colour_identity() +
    geom_point(color = "#FFFFB3") +
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
    ggtitle("Cluster 3") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc3) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 3 - females ############

# summary statistics
micro3_female %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro3_female %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro3_female %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro3_female %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc3 <- micro3_female %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc3 <- pwc3 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 3_female_diagnosis.pdf"))
ggplot(micro3_female, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#FFFFB3") + scale_colour_identity() +
    geom_point(color = "#FFFFB3") +
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
    ggtitle("Cluster 3 - Females") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc3) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 3 - males ############

# summary statistics
micro3_male %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro3_male %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro3_male %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro3_male %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc3 <- micro3_male %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc3 <- pwc3 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 3_male_diagnosis.pdf"))
ggplot(micro3_male, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#FFFFB3") + scale_colour_identity() +
    geom_point(color = "#FFFFB3") +
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
    ggtitle("Cluster 3 - Males") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc3) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 4 ############

# summary statistics
micro4 %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro4 %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro4 %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro4 %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc4 <- micro4 %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc4 <- pwc4 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 4_diagnosis.pdf"))
ggplot(micro4, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#80B1D3") + scale_colour_identity() +
    geom_point(color = "#80B1D3") +
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
    ggtitle("Cluster 4") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc4) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 4 - females ############

# summary statistics
micro4_female %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro4_female %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro4_female %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro4_female %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc4 <- micro4_female %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc4 <- pwc4 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 4_female_diagnosis.pdf"))
ggplot(micro4_female, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#80B1D3") + scale_colour_identity() +
    geom_point(color = "#80B1D3") +
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
    ggtitle("Cluster 4 - Females") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc4) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 4 - males ############

# summary statistics
micro4_male %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro4_male %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro4_male %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro4_male %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc4 <- micro4_male %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc4 <- pwc4 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 4_male_diagnosis.pdf"))
ggplot(micro4_male, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#80B1D3") + scale_colour_identity() +
    geom_point(color = "#80B1D3") +
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
    ggtitle("Cluster 4 - Males") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc4) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 5 ############

# summary statistics
micro5 %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro5 %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro5 %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro5 %>% anova_test(percentage ~ diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc5 <- micro5 %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc5 <- pwc5 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 5_diagnosis.pdf"))
ggplot(micro5, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#D9D9D9") + scale_colour_identity() +
    geom_point(color = "#D9D9D9") +
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
    ggtitle("Cluster 5") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc5) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 5 - females ############

# summary statistics
micro5_female %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro5_female %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro5_female %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro5_female %>% anova_test(percentage ~ diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc5 <- micro5_female %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc5 <- pwc5 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 5_female_diagnosis.pdf"))
ggplot(micro5_female, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#D9D9D9") + scale_colour_identity() +
    geom_point(color = "#D9D9D9") +
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
    ggtitle("Cluster 5 - Females") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc5) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 5 ############

# summary statistics
micro5_male %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro5_male %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro5_male %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro5_male %>% anova_test(percentage ~ diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc5 <- micro5_male %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc5 <- pwc5 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 5_male_diagnosis.pdf"))
ggplot(micro5_male, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#D9D9D9") + scale_colour_identity() +
    geom_point(color = "#D9D9D9") +
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
    ggtitle("Cluster 5 - Males") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc5) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()


######### Micro 6 ############

# summary statistics
micro6 %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro6 %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro6 %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro6 %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc6 <- micro6 %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc6 <- pwc6 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 6_diagnosis.pdf"))
ggplot(micro6, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#BC80BD") + scale_colour_identity() +
    geom_point(color = "#BC80BD") +
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
    ggtitle("Cluster 6") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc6) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 6 - females ############

# summary statistics
micro6_female %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro6_female %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro6_female %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro6_female %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc6 <- micro6_female %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc6 <- pwc6 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 6_female_diagnosis.pdf"))
ggplot(micro6, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#BC80BD") + scale_colour_identity() +
    geom_point(color = "#BC80BD") +
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
    ggtitle("Cluster 6 - Females") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc6) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 6 - males ############

# summary statistics
micro6_male %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro6_male %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro6_male %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro6_male %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc6 <- micro6_male %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc6 <- pwc6 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 6_male_diagnosis.pdf"))
ggplot(micro6, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#BC80BD") + scale_colour_identity() +
    geom_point(color = "#BC80BD") +
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
    ggtitle("Cluster 6 - Males") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc6) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 7 ############

# summary statistics
micro7 %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro7 %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro7 %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro7 %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc7 <- micro7 %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc7 <- pwc7 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 7_diagnosis.pdf"))
ggplot(micro7, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#FCCDE5") + scale_colour_identity() +
    geom_point(color = "#FCCDE5") +
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
    ggtitle("Cluster 7") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc7) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 7 ############

# summary statistics
micro7_female %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro7_female %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro7_female %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro7_female %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc7 <- micro7_female %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc7 <- pwc7 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 7_female_diagnosis.pdf"))
ggplot(micro7_female, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#FCCDE5") + scale_colour_identity() +
    geom_point(color = "#FCCDE5") +
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
    ggtitle("Cluster 7 - Females") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc7) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 7 - males ############

# summary statistics
micro7_male %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro7_male %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro7_male %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro7_male %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc7 <- micro7_male %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc7 <- pwc7 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 7_male_diagnosis.pdf"))
ggplot(micro7_male, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#FCCDE5") + scale_colour_identity() +
    geom_point(color = "#FCCDE5") +
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
    ggtitle("Cluster 7 - Males") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc7) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 8 ############

# summary statistics
micro8 %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro8 %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro8 %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro8 %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc8 <- micro8 %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc8 <- pwc8 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 8_diagnosis.pdf"))
ggplot(micro8, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#BEBADA") + scale_colour_identity() +
    geom_point(color = "#BEBADA") +
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
    ggtitle("Cluster 8") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc8) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 8 - females ############

# summary statistics
micro8_female %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro8_female %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro8_female %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro8_female %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc8 <- micro8_female %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc8 <- pwc8 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 8_female_diagnosis.pdf"))
ggplot(micro8_female, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#BEBADA") + scale_colour_identity() +
    geom_point(color = "#BEBADA") +
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
    ggtitle("Cluster 8 - Females") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc8) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 8 - Males ############

# summary statistics
micro8_male %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro8_male %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro8_male %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro8_male %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc8 <- micro8_male %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc8 <- pwc8 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 8_male_diagnosis.pdf"))
ggplot(micro8_male, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#BEBADA") + scale_colour_identity() +
    geom_point(color = "#BEBADA") +
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
    ggtitle("Cluster 8 - Males") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc8) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()


######### Micro 9 ############

# summary statistics
micro9 %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro9 %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro9 %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro9 %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc9 <- micro9 %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc9 <- pwc9 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 9_diagnosis.pdf"))
ggplot(micro9, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#FDB462") + scale_colour_identity() +
    geom_point(color = "#FDB462") +
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
    ggtitle("Cluster 9") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc9) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 9 - females ############

# summary statistics
micro9_female %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro9_female %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro9_female %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro9_female %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc9 <- micro9_female %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc9 <- pwc9 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 9_female_diagnosis.pdf"))
ggplot(micro9_female, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#FDB462") + scale_colour_identity() +
    geom_point(color = "#FDB462") +
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
    ggtitle("Cluster 9 - Females") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc9) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 9 - males ############

# summary statistics
micro9_male %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro9_male %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro9_male %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro9_male %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc9 <- micro9_male %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc9 <- pwc9 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 9_male_diagnosis.pdf"))
ggplot(micro9_male, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#FDB462") + scale_colour_identity() +
    geom_point(color = "#FDB462") +
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
    ggtitle("Cluster 9 - Males") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc9) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()


######### Micro 10 ############

# summary statistics
micro10 %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro10 %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro10 %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro10 %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc10 <- micro10 %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc10 <- pwc10 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 10_diagnosis.pdf"))
ggplot(micro10, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#FB8072") + scale_colour_identity() +
    geom_point(color = "#FB8072") +
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
    ggtitle("Cluster 10") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc10) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 10 - females ############

# summary statistics
micro10_female %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro10_female %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro10_female %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro10_female %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc10 <- micro10_female %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc10 <- pwc10 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 10_female_diagnosis.pdf"))
ggplot(micro10_female, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#FB8072") + scale_colour_identity() +
    geom_point(color = "#FB8072") +
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
    ggtitle("Cluster 10 - Females") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc10) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()


######### Micro 10 - Males ############

# summary statistics
micro10_male %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro10_male %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro10_male %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro10_male %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc10 <- micro10_male %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc10 <- pwc10 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 10_male_diagnosis.pdf"))
ggplot(micro10_male, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#FB8072") + scale_colour_identity() +
    geom_point(color = "#FB8072") +
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
    ggtitle("Cluster 10 - Males") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc10) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 11 ############

# summary statistics
micro11 %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro11 %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro11 %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro11 %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc11 <- micro11 %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc11 <- pwc11 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 11_diagnosis.pdf"))
ggplot(micro11, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#CCEBC5") + scale_colour_identity() +
    geom_point(color = "#CCEBC5") +
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
    ggtitle("Cluster 11") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc11) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 11 - females ############

# summary statistics
micro11_female %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro11_female %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro11_female %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro11_female %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc11 <- micro11_female %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc11 <- pwc11 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 11_female_diagnosis.pdf"))
ggplot(micro11_female, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#CCEBC5") + scale_colour_identity() +
    geom_point(color = "#CCEBC5") +
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
    ggtitle("Cluster 11 - Females") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc11) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 11 - males ############

# summary statistics
micro11_male %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro11_male %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro11_male %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro11_male %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc11 <- micro11_male %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc11 <- pwc11 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 11_male_diagnosis.pdf"))
ggplot(micro11_male, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#CCEBC5") + scale_colour_identity() +
    geom_point(color = "#CCEBC5") +
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
    ggtitle("Cluster 11 - Males") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc11) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 12 ############

# summary statistics
micro12 %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro12 %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro12 %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro12 %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc12 <- micro12 %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc12 <- pwc12 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 12_diagnosis.pdf"))
ggplot(micro12, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#FFED6F") + scale_colour_identity() +
    geom_point(color = "#FFED6F") +
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
    ggtitle("Cluster 12") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc12) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 12 - females ############

# summary statistics
micro12_female %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro12_female %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro12_female %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro12_female %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc12 <- micro12_female %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc12 <- pwc12 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 12_female_diagnosis.pdf"))
ggplot(micro12_female, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#FFED6F") + scale_colour_identity() +
    geom_point(color = "#FFED6F") +
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
    ggtitle("Cluster 12 - Females") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc12) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()

######### Micro 12 - males ############

# summary statistics
micro12_male %>% group_by(diagnosis) %>%
    get_summary_stats(percentage, type = "mean_sd")
# check outliers
micro12_male %>%
    group_by(diagnosis) %>%
    identify_outliers(percentage)
# normality assumption
  micro12_male %>%
    group_by(diagnosis) %>%
    shapiro_test(percentage) # p moet groter dan 0.05 zijn dan is het normaal verdeeld
# anova calculation
  res.aov <- micro12_male %>% anova_test(percentage ~  diagnosis)
  get_anova_table(res.aov) # ges is the amount of variablity due to the within-subjects factor

# post hoc pairwise comparisons
  pwc12 <- micro12_male %>% tukey_hsd(percentage ~  diagnosis)

# plot results
  # Visualization: box plots with p-values
  pwc12 <- pwc12 %>% add_xy_position(x = "percentage")

lineWidth <- 1
pointSize <- 20

pdf(file = paste0(Output.dir.subset.results,"Cluster distributions - Micro 12_male_diagnosis.pdf"))
ggplot(micro12_male, aes(x=diagnosis, y=percentage)) +
    geom_boxplot(color = "#FFED6F") + scale_colour_identity() +
    geom_point(color = "#FFED6F") +
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
    ggtitle("Cluster 12 - Males") + 
    xlab("Diagnosis") +
    ylab("Proportion of nuclei in cluster (%)") +
stat_pvalue_manual(pwc12) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE)
    )
dev.off()





