######--------Human pathway analysis--------######
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(limma)
library(openxlsx)
library(dplyr)
library(tidyverse)
library(data.table)
library(stringr)
library(fgsea)
library(enrichplot)
library(pheatmap)
library(dendextend)

#Load data for enrichment analysis
file_prefix <- 'xxx'
file_export <- 'xxx'
pathway_data <- read.xlsx(paste0('results/',file_prefix, 'DEGene_statistics_pval05.xlsx'))

#Filter genes 
DEgenes <- pathway_data %>% filter(pvalue < 0.05)  %>% pull('gene')
#DEgenes <- mapIds(org.Hs.eg.db, DEgenes, 'ENTREZID', 'SYMBOL')
DEG_Log2FC <-  pathway_data %>% filter(pvalue < 0.05) %>%  pull('log2FoldChange')
names(DEG_Log2FC) = DEgenes
DEG_Log2FC <- sort(DEG_Log2FC, decreasing = TRUE)
DEG_Log2FC

# DEG_Log2FC<-na.omit(DEG_Log2FC)
# DEG_Log2FC<-DEG_Log2FC[!is.na(DEG_Log2FC)]
# DEG_Log2FC

############ --------- FGSEA analysis ####
#Fgsea plot theme
theme_gsea <- theme(
  text = element_text(size = pointSize, colour = "black"),
  rect = element_blank(),
  line = element_line(size = lineWidth, colour = "black"),
  plot.title  = element_text(color="black", size = pointSize),
  axis.title  = element_text(size = pointSize * 0.8, colour = "black"),
  axis.text.x  = element_text(size = pointSize * 0.8, colour = "black"),
  axis.text.y  = element_text(size = pointSize * 0.8, colour = "black"),
  axis.ticks.x = element_line(size = 0.5, colour = "black"),
  axis.ticks.y = element_line(size = 0.5, colour = "black"),
  axis.line = element_line(size = lineWidth, colour = "black"),
  axis.line.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.title = element_text(size = pointSize * 0.7, colour = "black"),
  legend.text = element_text(size = pointSize * 0.4, colour = "black"),
  legend.key.height = unit(0.2, "cm"),
  legend.key.width = unit(0.2, "cm"))

#### CURATED pathways #####
pathways.all.curated <- gmtPathways("GSEA_database/c2.all.v2022.1.Hs.entrez.gmt") 

fgseaRes <- fgsea(pathways=pathways.all.curated, stats=DEG_Log2FC)

fgseaResTidy_curated <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  filter(pval<0.05)

write.xlsx(fgseaResTidy_curated, paste0('results/Pathway_analysis/FGSEA_enrichment/',file_export, 'Curated_pathways.xlsx'))

##### FILTER SPECIFIC PATHWAYS ####
fgseaResTidy_curated_pval05 <- fgseaResTidy_curated %>% filter(pval < 0.1)
fgseaResTidy_curated_pval05 <- fgseaResTidy_curated_pval05 %>% filter(pathway %in% c('REACTOME_NEUTROPHIL_DEGRANULATION','LIAN_NEUTROPHIL_GRANULE_CONSTITUENTS','MARTINELLI_IMMATURE_NEUTROPHIL_UP',
                                                                                     'MAHAJAN_RESPONSE_TO_IL1A_UP', 'REACTOME_FCGR3A_MEDIATED_IL10_SYNTHESIS', 'KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY',
                                                                                     'BIOCARTA_IL17_PATHWAY','WP_IL17_SIGNALING_PATHWAY','MAHAJAN_RESPONSE_TO_IL1A_UP','REACTOME_COLLAGEN_DEGRADATION',
                                                                                     'COULOUARN_TEMPORAL_TGFB1_SIGNATURE_UP','WP_ALLOGRAFT_REJECTION','REACTOME_SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX',
                                                                                     'VERRECCHIA_RESPONSE_TO_TGFB1_C2', 'REACTOME_REGULATION_OF_IFNA_IFNB_SIGNALING','MARZEC_IL2_SIGNALING_UP',
                                                                                     'GILMORE_CORE_NFKB_PATHWAY','KEGG_LYSOSOME','BIOCARTA_NEUTROPHIL_PATHWAY','REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS',
                                                                                     'WP_DEGRADATION_PATHWAY_OF_SPHINGOLIPIDS_INCLUDING_DISEASES','LIU_IL13_PRIMING_MODEL',
                                                                                     'WP_IL4_SIGNALING_PATHWAY', 'REACTOME_SIGNALING_BY_TGFB_FAMILY_MEMBERS','REACTOME_SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX',
                                                                                     'KONDO_HYPOXIA', 'REACTOME_REGULATED_NECROSIS','REACTOME_INTERLEUKIN_4_AND_INTERLEUKIN_13_SIGNALING',
                                                                                     'HARRIS_HYPOXIA','REACTOME_NEUTROPHIL_DEGRANULATION','KEGG_INSULIN_SIGNALING_PATHWAY','BIOCARTA_NFAT_PATHWAY'))
fgseaResTidy_curated_pval05$log_pval <- -log10(fgseaResTidy_curated_pval05$pval)
barplot_curated <- ggplot(fgseaResTidy_curated_pval05, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES), color = 'black') +
  coord_flip() +
  theme_gsea +
  scale_fill_gradient2(low = 'blue', mid = "white",
                       high = 'red', midpoint = 0, space = "rgb",
                       na.value = "grey50", guide = "colourbar") +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Curated pathways NES from GSEA") + 
  geom_vline(xintercept = 0, color ='black')
barplot_curated

ggsave(paste0('plots/pathways/FGSEA_pathways/', file_export,'curated_pathways.pdf'),
       barplot_curated,
       dpi = 300)

selected_curated_pathways <-  c('REACTOME_NEUTROPHIL_DEGRANULATION','LIAN_NEUTROPHIL_GRANULE_CONSTITUENTS','MARTINELLI_IMMATURE_NEUTROPHIL_UP',
                                'MAHAJAN_RESPONSE_TO_IL1A_UP', 'REACTOME_FCGR3A_MEDIATED_IL10_SYNTHESIS', 'KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY',
                                'BIOCARTA_IL17_PATHWAY','WP_IL17_SIGNALING_PATHWAY','MAHAJAN_RESPONSE_TO_IL1A_UP','REACTOME_COLLAGEN_DEGRADATION',
                                'COULOUARN_TEMPORAL_TGFB1_SIGNATURE_UP','WP_ALLOGRAFT_REJECTION','REACTOME_SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX',
                                'VERRECCHIA_RESPONSE_TO_TGFB1_C2', 'REACTOME_REGULATION_OF_IFNA_IFNB_SIGNALING','MARZEC_IL2_SIGNALING_UP',
                                'GILMORE_CORE_NFKB_PATHWAY','KEGG_LYSOSOME','BIOCARTA_NEUTROPHIL_PATHWAY','REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS',
                                'WP_DEGRADATION_PATHWAY_OF_SPHINGOLIPIDS_INCLUDING_DISEASES','LIU_IL13_PRIMING_MODEL',
                                'WP_IL4_SIGNALING_PATHWAY', 'REACTOME_SIGNALING_BY_TGFB_FAMILY_MEMBERS','REACTOME_SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX',
                                'KONDO_HYPOXIA', 'REACTOME_REGULATED_NECROSIS')
for (i in selected_curated_pathways) {
  print(i)
  genes_curated <- pathways.all.curated[[i]]
  genes_curated_mapped <- bitr(genes_curated, fromType = "ENTREZID",
                               toType = "SYMBOL",
                               OrgDb = org.Hs.eg.db)
  write.xlsx(genes_curated_mapped, paste0('results/Pathway_analysis/FGSEA_enrichment/pathway_genes/',i,'_','curated_pathways.xlsx'))
}


#Custom Neutrophil pathway
pathways.custom <- gmtPathways("GSEA_database/custom_neutrophil_pathways/neutrophil_custom_singlecell_test.v2023.1.Hs.gmt") 
fgseaRes <- fgsea(pathways=pathways.custom, stats=DEG_Log2FC)
fgseaResTid_custom <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  filter(pval<0.1)

write.xlsx(fgseaResTid_custom, paste0('results/Pathway_analysis/FGSEA_enrichment/',file_export, '_neutro_singlecell_custom_pathways.xlsx'))


pathways.custom <- gmtPathways("GSEA_database/custom_neutrophil_pathways/neutrophil_custom_singlecell.v2023.1.Hs.gmt") 
fgseaRes <- fgsea(pathways=pathways.custom, stats=DEG_Log2FC)

fgseaResTid_custom <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  filter(pval<0.1)

write.xlsx(fgseaResTid_custom, paste0('results/Pathway_analysis/FGSEA_enrichment/',file_export, '_neutro_singlecell_custom_pathways.xlsx'))

#Plot enrichment - Set up function 
lineColor <- 'firebrick'
plotCustomEnrichment <- function(pathway, stats,
                                 gseaParam=1,
                                 ticksSize=0.2) {
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
  statsAdj <- statsAdj / max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))
  diff <- (max(tops) - min(bottoms)) / 8
  # Getting rid of NOTEs
  x=y=NULL
  g <- ggplot(toPlot, aes(x=x, y=y)) +
    #geom_point(color='black', size=4) +
    geom_point(fill=lineColor, size=3, pch = 21) +
    #geom_hline(yintercept=max(tops), colour="red", linetype="dashed") +
    #geom_hline(yintercept=min(bottoms), colour="red", linetype="dashed") +
    geom_hline(yintercept=0, colour="black") +
    #geom_line(color=lineColor, size =1) + 
    theme_bw() +
    geom_segment(data=data.frame(x=pathway),
                 mapping=aes(x=x, y=-diff/2,xend=x, yend=diff/2),size=ticksSize) +
    theme(axis.line = element_line(size = 0.5, colour = "black"),
          axis.ticks.x = element_line(size = 0.5, colour = "black"),
          axis.ticks.y = element_line(size = 0.5, colour = "black"),
          axis.text.x  = element_text(size = 20, colour = "black"),
          axis.text.y  = element_text(size = 20, colour = "black"),
          plot.title  = element_text(color="black", size = 30),
          axis.title  = element_text(size = 20, colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    labs(x="rank", y="enrichment score")
  g
}

Enrichment_plot <- plotCustomEnrichment(pathways.custom[["Mature_Neutrophils"]], stats=DEG_Log2FC) + labs(title="Mature Neutrophils - MCI/AD vs Ctrl") + labs(subtitle=paste0('NES =',fgseaResTid_custom$NES,"   ",'P-adjusted =', fgseaResTid_custom$padj))
Enrichment_plot
ggsave(paste0('plots/pathways/FGSEA_pathways/Custom_single_cell_neutro_patway_color_', file_export,'.pdf'),
       Enrichment_plot,
       dpi = 300)



###### AGING RELATED SIGNATURE NEUTROPHILS #######ö
file_prefix <- 'xxx'
file_export <- 'xxx'

pathway_data <- read.xlsx(paste0('results/',file_prefix, 'DEGene_statistics_pval05.xlsx'))
DEgenes <- pathway_data %>% filter(pvalue < 0.05)  %>% pull('gene')
DEG_Log2FC <-  pathway_data %>% filter(pvalue < 0.05) %>%  pull('log2FoldChange')

names(DEG_Log2FC) = DEgenes
DEG_Log2FC <- sort(DEG_Log2FC, decreasing = TRUE)
DEG_Log2FC

#GSEA Aging associated gene pathway analysis
pathways.aging_up <- gmtPathways("Aging_UP_genesets.v2023.1.Hs.gmt") 
fgseaRes <- fgsea(pathways=pathways.aging_up, stats=DEG_Log2FC)

fgseaResTidy_aging <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  filter(pval<0.05)

write.xlsx(fgseaResTidy_aging, paste0('results/Pathway_analysis/FGSEA_enrichment/',file_export, 'aging_pathways.xlsx'))


#Extract aging related genes and make heatmap
toMatch <- c("AGING")
aging_data <- dplyr::filter(fgseaResTidy_aging,pathway %in% unique(grep(paste(toMatch,collapse="|"), 
                                                                        fgseaResTidy_aging$pathway, value=TRUE))) 
#Filter aging related pathways
aging_data <- as.list(aging_data$leadingEdge)
aging_data_list <- aging_data[[1]]
aging_data_list

#Create Heatmap for aging data
aging_heatmap_data <- read.xlsx(paste0('results/',aging_pathway_data_prefix, 'DEGene_counts_pval05.xlsx'))
aging_heatmap_data <- aging_heatmap_data %>% filter(gene %in% aging_data_list) %>% column_to_rownames('gene')

# Define a normalization function to calculate Z scores.
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)}
# Normalize your data according to Z score using cal_z_score function.
aging_heatmap_data_norm <- t(apply(aging_heatmap_data, 1, cal_z_score))
#### HEATMAP  ####
cluster_no = 3
gaps = c(11)
# Create heatmap
pdf(file = paste0('plots/pathways/FGSEA_pathways/', 'Aging_genes_',file_prefix, 'pval05.pdf'), pointsize = 10)
phm_full <- pheatmap(aging_heatmap_data_norm,
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = F,
                     cutree_row = cluster_no,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 10,
                     treeheight_col = 50,
                     treeheight_row = 50,
                     cellwidth = 10,
                     scale = 'row',
                     cellheight = 10
)
dev.off() 

