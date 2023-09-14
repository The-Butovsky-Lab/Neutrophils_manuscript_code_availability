#Loading packagees
library(tidyverse)
library(pheatmap)
library(dendextend)
library(DESeq2)
library(openxlsx)
library(dplyr)

# Define a normalization function to calculate Z scores.
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)}

#Analyzis parameters
file_prefix = 'xxx'

# Import in data and subset it for how many genes you want to see.
data <- read.xlsx(paste0('results/', file_prefix, 'DEGene_counts_pval05.xlsx'), rowNames = T) 

#### SET UP DATA ####
data_final <- data %>%
  rownames_to_column('gene') %>%
  #filter(!str_detect(gene, c(unwanted_genes))) %>%
  #anti_join(unwanted_clusters, by = 'gene') %>%
  column_to_rownames('gene')
data_subset_200 <- data_final %>% head(n = 200)
data_subset_400 <- data_final %>% head(n = 400)

data_grouped <- data.frame(
  xxx = rowMeans(dplyr::select(data, contains(("xxx")))))
data_grouped_subset_400 <- data_grouped %>% head(n = 400)

# Normalize your data according to Z score using cal_z_score function.
data_norm <- t(apply(data_final, 1, cal_z_score))
data_subset_200 <- t(apply(data_subset_200, 1, cal_z_score))
data_subset_400 <- t(apply(data_subset_400, 1, cal_z_score))
data_grouped_norm <- t(apply(data_grouped, 1, cal_z_score))
data_grouped_subset_400_norm <- t(apply(data_grouped_subset_400, 1, cal_z_score))


#### HEATMAP  ####
cluster_no = 6
gaps = c(5,11,18)
# Create heatmap
pdf(file = paste0('plots/heatmaps/', file_prefix, 'pval05.pdf'), pointsize = 10, width = 15, height = 40)
phm_full <- pheatmap(data_norm,
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = cluster_no,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = F,
                     border_color = 'NA',
                     fontsize = 10,
                     treeheight_col = 50,
                     treeheight_row = 50,
                     #cellwidth = 25,
                     scale = 'row'
                     #cellheight = 3
                     )
dev.off() 


pdf(file = paste0('plots/heatmaps/', file_prefix, 'labelled_full_pval05.pdf'), pointsize = 10, width = 15, height = 180)
phm_full <- pheatmap(data_norm,
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = cluster_no,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = NA,
                     fontsize = 6,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     #cellwidth = 25,
                     scale = 'row',
                     cellheight = 7)
dev.off() 


pdf(file = paste0('plots/heatmaps/', file_prefix, '200.pdf'), pointsize = 10, height = 32, width = 15)
phm_100 <- pheatmap(data_subset_200,
                    color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                    breaks = seq(-2, 2, by = 0.1),
                    kmeans_k = NA,
                    cluster_rows = T,
                    cutree_row = cluster_no,
                    cluster_cols = F,
                    #cutree_cols = 4,
                    gaps_col = gaps,              
                    legend = TRUE,
                    show_rownames = T,
                    #cellwidth = 25,
                    cellheight = 9,
                    treeheight_col = 0,
                    treeheight_row = 50,
                    border_color = 'NA',
                    fontsize = 8,
                    scale = 'row')
dev.off()



pdf(file = paste0('plots/heatmaps/', file_prefix, 'heatmap_400.pdf'), pointsize = 10, height = 65, width = 15)
phm_full <- pheatmap(data_subset_400,
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = cluster_no,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     #cellwidth = 25,
                     scale = 'row',
                     cellheight = 9)
dev.off()

#Padj
# pdf(file = paste0('plots/', file_prefix, 'padj05.pdf'), pointsize = 10, height = 32, width = 15)
# phm_100 <- pheatmap(data_norm,
#                     color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
#                     breaks = seq(-2, 2, by = 0.1),
#                     kmeans_k = NA,
#                     cluster_rows = T,
#                     cutree_row = 10,
#                     cluster_cols = F,
#                     #cutree_cols = 4,
#                     gaps_col = c(),              
#                     legend = TRUE,
#                     show_rownames = T,
#                     #cellwidth = 25,
#                     cellheight = 9,
#                     treeheight_col = 0,
#                     treeheight_row = 50,
#                     border_color = 'NA',
#                     fontsize = 8,
#                     scale = 'row'
# )
# dev.off()

### GOURPED
#Labelled Grouped
pdf(file = paste0('plots/heatmaps/', file_prefix, 'grouped_labelled_pval.pdf'), height = 120, width = 8)
phm_full <- pheatmap(data_grouped_norm,
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = cluster_no,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 5,
                     treeheight_col = 50,
                     treeheight_row = 50,
                     cellwidth = 10,
                     cellheight = 6,
                     scale = 'row')
dev.off()

#Grouped no labels
pdf(file = paste0('plots/heatmaps/', file_prefix, 'grouped_pval.pdf'), height = 90, width = 8)
phm_full_grouped <- pheatmap(data_grouped_norm,
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = cluster_no,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = F,
                     border_color = 'NA',
                     fontsize = 5,
                     treeheight_col = 50,
                     treeheight_row = 50,
                     #cellwidth = 10,
                     #cellheight = 6,
                     scale = 'row',
                     clustering_distance_rows = "euclidean")
dev.off()


pdf(file = paste0('plots/heatmaps/', file_prefix, 'grouped_400.pdf'), pointsize = 10, height = 70)
phm_full <- pheatmap(data_grouped_subset_400_norm,
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = cluster_no,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 50,
                     treeheight_row = 50,
                     #cellwidth = 25,
                     scale = 'row',
                     cellheight = 9
)
dev.off()

#### CLUSTER EXTRACTION ####
# Build cluster tables from dendrogram based on k
gene_clusters <- data.frame(sort(cutree(phm_full_grouped$tree_row, k=6))) %>%
  rownames_to_column(var="gene") 

temp = gene_clusters$sort.cutree.phm_full_grouped.tree_row..k...6..

gene_clusters <- gene_clusters %>% 
  dplyr::mutate(cluster = temp) %>%
  dplyr::select(gene, cluster)

write.xlsx(gene_clusters, file=paste0('results/', file_prefix, 'grouped_6clusters_clustergenes.xlsx'))

#### HEATMAP CUSTOMIZATION #####
# Subset and organize clusters (iterative process)
gene_clusters_clean <- gene_clusters %>%
  arrange(match(cluster, c(1,6,8,3,4,7,5,2)))


# Reorder transcript count data to match custom order
#data_ordered <- data_final[gene_clusters_clean$gene, ]    
data_ordered <- data_final[gene_clusters_clean$gene,]

# Build annotation dataframes metadata
gene_clusters_clean <- gene_clusters_clean %>%
  mutate(primary_effect = ifelse(cluster == 3, 'xx',
                            ifelse(cluster == 7, 'xx',
                            ifelse(cluster %in% c(1,8),'xx',
                            ifelse(cluster == 4, 'xx',
                            ifelse(cluster == 5, 'xx',
                            ifelse(cluster %in% c(2), 'xx',
                                 'xx')))))))


# Subset out annotation data for plotting
cluster_anno <- gene_clusters_clean %>% select(primary_effect)
rownames(cluster_anno) <- rownames(data_ordered)

##
data_ordered_grouped <- data.frame(
  xxx = rowMeans(select(data_ordered, contains(("xxx")))))


#Set up color scheme
cols <- brewer.pal(6, "Paired")
cols
cols_heatmap_annotation <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C")

# Create heatmap
pdf(file = paste0('plots/heatmaps/', file_prefix, 'custom_grouped_pval05.pdf'), pointsize = 10)
phm_custom <- pheatmap(data_ordered_grouped,
                       color = colorRampPalette(c("navyblue", "white", "firebrick3"))(41),
                       breaks = seq(-2, 2, by = 0.1),
                       kmeans_k = NA,
                       cluster_rows = F,
                       #cutree_rows = 3,
                       cluster_cols = F, 
                       gaps_row = c(172,364,443,796,975,1045), #Male
                       #gaps_row = c(142,216,395,514,600,755), #Female 
                       #gaps_row = c(142,216,386,395,514,600,755),    
                       gaps_col = c(),     
                       annotation_legend = T,
                       annotation_row = cluster_anno,
                       annotation_names_row = F,
                       #annotation_col = column_anno,
                       annotation_names_col = F,
                       #annotation_colors = cols_heatmap_annotation,
                       show_rownames = F,
                       show_colnames = T,
                       border_color = NA,
                       cellwidth = 25,
                       treeheight_col = 0,
                       treeheight_row = 0,
                       fontsize_row = 7,
                       scale = 'row',
                       legend = TRUE
)
dev.off()

#Custom labelled
pdf(file = paste0('plots/heatmaps/', file_prefix, 'custom_grouped_labelled_pval05.pdf'),height =150, width = 50)
phm_custom <- pheatmap(data_ordered_grouped,
                       color = colorRampPalette(c("navyblue", "white", "firebrick3"))(41),
                       breaks = seq(-2, 2, by = 0.1),
                       kmeans_k = NA,
                       cluster_rows = F,
                       #cutree_rows = 3,
                       cluster_cols = F,
                       gaps_row = c(172,364,443,796,975,1045), #Male
                       #gaps_row = c(142,216,395,514,600,755),
                       #gaps_row = c(142,216,386,395,514,600,755),    
                       gaps_col = c(),     
                       annotation_legend = T,
                       annotation_row = cluster_anno,
                       annotation_names_row = F,
                       #annotation_col = column_anno,
                       annotation_names_col = F,
                       #annotation_colors = anno_colors,
                       show_rownames = T,
                       show_colnames = T,
                       border_color = NA,
                       cellwidth = 25,
                       treeheight_col = 0,
                       treeheight_row = 0,
                       fontsize_row = 8,
                       scale = 'row',
                       legend = TRUE,
                       cellheight = 9
)
dev.off()
