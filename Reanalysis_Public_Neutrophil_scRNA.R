if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")

list.of.packages <- c('tidyverse', 'Seurat', 'Signac', 'enrichR', 'openxlsx', 'patchwork', 
                      'data.table', 'dittoSeq', 'ggplot2','RColorBrewer','limma',
                      'openxlsx','dplyr','data.table','stringr','fgsea', 'enrichplot',
                      'org.Hs.eg.db', 'clustree')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

lapply(list.of.packages, require, character.only = TRUE)
options(future.globals.maxSize = 4000 * 1024^5)

###### Re-annalysis SCHUTLE-SCHREPPING et al. 2020 ######
seuratwb <- readRDS('cibersort/single_cell/Schulte_Schrepping/seurat_objects/seurat_COVID19_freshWB-PBMC_cohort2_rhapsody_jonas_FG_2020-08-18.rds')
barcodecluster <- as.character(seuratwb$RNA_snn_res.0.8)
temp <- seuratwb$RNA_snn_res.0.8
DimPlot(temp)
names(barcodecluster) <- names(temp)
barcodecluster <- cbind(rownames(barcodecluster),barcodecluster)
barcodecluster <- as.data.frame(barcodecluster)
barcodecluster$barcode <- as.character(rownames(barcodecluster))
barcodecluster$barcode[barcodecluster$barcodecluster %in% c("0","1","4","8")] <- "Mature_Neutrophil"
barcodecluster$barcode[barcodecluster$barcodecluster %in% c("12","15")] <- "Immature_Neutrophil"
barcodecluster$barcode[barcodecluster$barcodecluster %in% c("3","10","14","18")] <- "Monocyte"
barcodecluster$barcode[barcodecluster$barcodecluster %in% c("2","5","6","13","16")] <- "T_NK"
barcodecluster$barcode[barcodecluster$barcodecluster %in% c("7","22")] <- "B"
barcodecluster$barcode[barcodecluster$barcodecluster %in% c("19")] <- "Plasmablast"
barcodecluster$barcode[barcodecluster$barcodecluster %in% c("9","11","17","20","21","23","24")] <- "Other"
seuratwb <- AddMetaData(seuratwb, metadata = subset(barcodecluster, select = c("barcode")), col.name = "Cell.Type")


Idents(seuratwb) <- 'Cell.Type'
levels(seuratwb)
DimPlot(seuratwb)
seuratwb_exprTable <- AverageExpression(seuratwb)
seuratwb_exprTable <- seuratwb_exprTable[['RNA']] %>% as.data.frame() 

##### Subsetting for only neutrophils #####
Idents(seuratwb) <- 'RNA_snn_res.0.8'
levels(seuratwb)
DimPlot(seuratwb, label = T)
dev.off()
Neutrophils <- subset(seuratwb, idents = c(1,0,12,15,4))

#Recluster Neutrophils
DimPlot(Neutrophils)
Idents(Neutrophils) <- 'RNA_snn_res.0.8'
DimPlot(Neutrophils, label = T)

################################## RECLUSTER FILTERED NEUTROPHILS #######################################
#RECLUSTER FOR NEUTROPHILS
dataset_filtered_reclust <- ScaleData(Neutrophils)
dataset_filtered_reclust <- FindVariableFeatures(dataset_filtered_reclust, selection.method = 'vst')
dataset_filtered_reclust <- RunPCA(dataset_filtered_reclust, features = VariableFeatures(object = dataset_filtered_reclust))
ElbowPlot(dataset_filtered_reclust)
# Determine percent of variation associated with each PC
pct <- dataset_filtered_reclust[["pca"]]@stdev / sum(dataset_filtered_reclust[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
pcs
# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))
# Elbow plot to visualize 
elbow_rank_vs_pcs <- ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
pdf(file = paste0('QC/dataset_filtered_reclust_elbow_rank_vs_pcs.pdf'), pointsize = 10, width = 16)
elbow_rank_vs_pcs
dev.off()

#Heatmap exploration
pdf(file = paste0('QC/dataset_filtered_reclust_heatmap_PCreductions_', projectname, '.pdf'), width = 20, height = 30)
DimHeatmap(dataset_filtered_reclust, dims = 1:20, cells = 1000, balanced = TRUE)
dev.off()

# Cluster cells according to elbow plot dimension choice
dataset_filtered_reclust <- FindNeighbors(dataset_filtered_reclust, dims = 1:9)

#Investigate resolution
dataset_filtered_reclust <- FindClusters(dataset_filtered_reclust, resolution = 1.2)
dataset_filtered_reclust <- FindClusters(dataset_filtered_reclust, resolution = 1.1)
dataset_filtered_reclust <- FindClusters(dataset_filtered_reclust, resolution = 1)
dataset_filtered_reclust <- FindClusters(dataset_filtered_reclust, resolution = 0.9)
dataset_filtered_reclust <- FindClusters(dataset_filtered_reclust, resolution = 0.8)
dataset_filtered_reclust <- FindClusters(dataset_filtered_reclust, resolution = 0.7)
dataset_filtered_reclust <- FindClusters(dataset_filtered_reclust, resolution = 0.6)
dataset_filtered_reclust <- FindClusters(dataset_filtered_reclust, resolution = 0.5)
dataset_filtered_reclust <- FindClusters(dataset_filtered_reclust, resolution = 0.4)
dataset_filtered_reclust <- FindClusters(dataset_filtered_reclust, resolution = 0.3)
dataset_filtered_reclust <- FindClusters(dataset_filtered_reclust, resolution = 0.2)

#Clustree exploration
cluster_tree <- clustree(dataset_filtered_reclust)
pdf(file = paste0('QC/dataset_filtered_reclust_cluster_tree_', projectname, '.pdf'), width = 12, height = 9)
cluster_tree
dev.off()

dataset_filtered_reclust <- FindNeighbors(dataset_filtered_reclust, dims = 1:9)
dataset_filtered_reclust <- FindClusters(dataset_filtered_reclust, resolution = 0.2)
dataset_filtered_reclust <- RunUMAP(dataset_filtered_reclust, dims = 1:9)

##Save seurat object
saveRDS(dataset_filtered_reclust,'Neutrophils_dim9_res02.rds')


#Find cluster markers for neutrophils
Cluster_markers_dataset_filtered_reclust <- FindAllMarkers(dataset_filtered_reclust, only.pos = T, min.pct = 0.1, logfc.threshold = 0.1)
TOP_Cluster_markers_dataset_filtered_reclust <- Cluster_markers_dataset_filtered_reclust %>% group_by(cluster) %>% top_n(n=100, wt = avg_log2FC)

# Cluster mapping
Idents(dataset_filtered_reclust) <- 'seurat_clusters'
levels(dataset_filtered_reclust)
new_mappings <- c('mature', 'Regulatory_INF_Type1', 'Intermediate', 'Immature')
names(new_mappings) <- levels(dataset_filtered_reclust)
dataset_filtered_reclust <- RenameIdents(dataset_filtered_reclust, new_mappings)
dataset_filtered_reclust$Neutrophil_states_4 <- Idents(dataset_filtered_reclust)
DimPlot(dataset_filtered_reclust, label = T)

#DimPlot
group_by_dimplot <- 'Neutrophil_states_4'
pdf(file = paste0('cibersort/single_cell/Schulte_Schrepping/plots/umap_dataset_filtered_reclust_neutrophils','_groupby_', group_by_dimplot,'.pdf'), width = 12, height = 12)
DimPlot(dataset_filtered_reclust, label = T, group.by = group_by_dimplot)
dev.off()

#no group DimPlot
pdf(file = paste0('cibersort/single_cell/Schulte_Schrepping/plots/umap_dataset_filtered_reclust_neutrophils_4subtypes_res02','.pdf'), width = 12, height = 12)
DimPlot(dataset_filtered_reclust, label = T)
dev.off()

#### DOTPLOT OF MARKER GENES #####
Dotplot_data <- DotPlot(dataset_filtered_reclust, features = c( 
  'IL18','CAMP','MMP8','LGALS3', 'LYZ', 'TSPO','GRN', 'ITGAM','MMP9','TGFBR1','S100A12', #Immature
  'S100A6','IL17RA', 'IL18R1', 'IL18RAP', 'IL4R',  'CR1', #Intermediate
  'IL10RB','TGFB1','CD44', 'CD37','CD177','ISG15', 'IFIT3', 'IFI44', 'SERPING1', 'IL1RN',  'IRF7',  'STAT1', 'LGALS9', 'CD274', 'IL2RG', 'NFKBIA',
  'CXCL8',  'CXCR4', 'TGFBR2', 'PTGS2', 'NR3C1','TREM1', 'FTH1' #Mature
  ))
Dotplot_data <-  as.data.frame(Dotplot_data$data)

dotplot_subtype_markers <- ggplot(Dotplot_data,aes(x = id, y = features.plot)) + 
  geom_point(aes(size = pct.exp),color ='black', stroke = 1.5) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled), stroke = 1) +
  theme(text = element_text(size = pointSize, colour = "black"),
        rect = element_blank(),
        line = element_line(size = lineWidth, colour = "black"),
        plot.title  = element_text(color="black", size=20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x  = element_text(size = pointSize , colour = "black", angle = 90, vjust=0.1),
        axis.text.y  = element_text(size = pointSize , colour = "black",face = 'italic' ),
        axis.ticks = element_line(size = lineWidth, colour = "black"),
        axis.line = element_line(size = lineWidth, colour = "black"),
        legend.position = "right",
        legend.title = element_text(size = pointSize , colour = "black"),
        legend.text = element_text(size = pointSize , colour = "black"),
        legend.key.height = unit(1, "cm"),
        legend.key.width = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(size = lineWidth, colour = "black"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0,size = pointSize),
        #strip.background = element_rect(size = 1))+
        strip.text.y = element_text(size=10,angle = 90)) +
  scale_color_gradient2(limits=c(-1.5,1.5), low="navyblue", mid="white", high = "firebrick3", na.value = 'red') +
  labs(title = str_wrap(paste0("Marker genes per Neutrophil Subcluster"), 40))
  #facet_grid(group~., scales = "free_y", space='free', switch = "both")
dotplot_subtype_markers

ggsave(
  paste0('cibersort/single_cell/Schulte_Schrepping/plots/', 'Dotplot_dataset_filtered_intermediate_4_SUBCELLTYES', '.pdf'),
  plot = dotplot_subtype_markers,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 170,
  height = 300,
  units = c( "mm"),
  dpi = 600
)


#### DOTPLOT OF FGSEA PATHWAY ACROSS SUBCLUSTER ####
#Canonical gene sets
pathway_data <- Cluster_markers_dataset_filtered_reclust
rm(canonical_data_total)
canonical_data_total <- data.frame()

pathways.canoical <- gmtPathways("cibersort/single_cell/Schulte_Schrepping/GSEA_pathways/c2.all.v2022.1.Hs.entrez.gmt") 

SubCellTypes_neutro <- pathway_data %>%  pull(cluster) %>% unique()
SubCellTypes_neutro

for (i in SubCellTypes_neutro){
  pathway_analysis <- pathway_data %>% filter(cluster == paste0(i))
  DEgenes <- pathway_analysis %>% filter(p_val < 0.05) %>% pull('gene')
  DEgenes <- mapIds(org.Hs.eg.db, DEgenes, 'ENTREZID', 'SYMBOL')
  DEG_Log2FC <-  pathway_analysis %>% filter(p_val < 0.05) %>%  pull('avg_log2FC')
  names(DEG_Log2FC) = DEgenes
  DEG_Log2FC <- sort(DEG_Log2FC, decreasing = TRUE)
  fgseaRes <- fgsea(pathways=pathways.canoical, stats=DEG_Log2FC, scoreType ='pos')
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES)) %>%
    filter(pval<0.05)
  fgseaResTidy$log_p <- -log10(fgseaResTidy$pval)
  fgseaResTidy$group <- paste0(i)
  canonical_data_total <- rbind(canonical_data_total, fgseaResTidy)
}
write.xlsx(canonical_data_total, paste0('cibersort/single_cell/Schulte_Schrepping/results/','Canonical_pathways.xlsx'))

#read selected pathways
canonical_selected_pathways <- read.xlsx(paste0('cibersort/single_cell/Schulte_Schrepping/results/Canonical_pathways.xlsx'), sheet = 2)
canonical_selected_pathways_filter <- canonical_selected_pathways %>% select(c('pathway','group','NES','pval'))

# Define a function to remove characters before the first underscore
remove_before_underscore <- function(string) {
  index <- regexpr("_", string)  
  if (index > 0) {
    substring(string, index + 1) 
  } else {
    string  
  }
}

# Apply the function to the 'strings' column
canonical_selected_pathways_filter$pathway <- sapply(canonical_selected_pathways_filter$pathway, remove_before_underscore)
canonical_selected_pathways_filter$pathway <- gsub("_CLUSTER_P6","",canonical_selected_pathways_filter$pathway)
canonical_selected_pathways_filter$pathway <- gsub("_CLUSTER_P4","",canonical_selected_pathways_filter$pathway)
canonical_selected_pathways_filter$pathway <- gsub("_4HR","",canonical_selected_pathways_filter$pathway)
canonical_selected_pathways_filter$pathway <- gsub("ALL_","",canonical_selected_pathways_filter$pathway)
canonical_selected_pathways_filter$pathway <- gsub("_"," ",canonical_selected_pathways_filter$pathway)

canonical_selected_pathways_filter$log_pval <- -log10(canonical_selected_pathways_filter$pval)

#Level pathway groups:
canonical_selected_pathways_filter$group <- factor(canonical_selected_pathways_filter$group, levels = c('mature','Intermediate','Regulatory_INF_Type1','Immature'))

#Makedotplot
lineWidth = 1
pointSize = 15
dotplot <- ggplot(canonical_selected_pathways_filter,aes(x = group, y = reorder(pathway,NES))) + 
  geom_point(aes(size = log_pval),color ='black', stroke = 4) +
  geom_point(aes(size = log_pval, color = NES), stroke = 3) +
  theme(text = element_text(size = pointSize, colour = "black"),
        rect = element_blank(),
        line = element_line(size = lineWidth, colour = "black"),
        plot.title  = element_text(color="black", size=20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x  = element_text(size = pointSize , colour = "black", angle = 90, vjust=0.1),
        axis.text.y  = element_text(size = pointSize , colour = "black"),
        axis.ticks = element_line(size = lineWidth, colour = "black"),
        axis.line = element_line(size = lineWidth, colour = "black"),
        legend.position = "right",
        legend.title = element_text(size = pointSize , colour = "black"),
        legend.text = element_text(size = pointSize , colour = "black"),
        legend.key.height = unit(1, "cm"),
        legend.key.width = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(size = lineWidth, colour = "black"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0,size = pointSize),
        strip.text.y = element_text(size=10,angle = 90)) +
  scale_color_gradientn(colours = c("firebrick2","firebrick", "red2", "red", "orange","purple"),
                        values = c(2,1.5,1,0.5,0)) +
  labs(title = str_wrap(paste0("FGSEA Pathways per Neutrophil Subcluster"), 40)) +
  facet_grid(group~., scales = "free_y", space='free', switch = "both")
dotplot

ggsave(paste0('cibersort/single_cell/Schulte_Schrepping/plots/pathways/FGSEA/Doplot_per_4subcelltypes_res02_','canonical_pathways_selected','.pdf'),
       dotplot,
       units = c("mm"),
       width = 300,
       height = 300,
       dpi = 300)














## Heatmap
canonical_selected_pathways <- read.xlsx(paste0('results/canonical_data_total.xlsx'), sheet = 2)%>% pull(pathway)
canonical_data_selected <- canonical_data_total %>% filter(pathway %in% canonical_selected_pathways)

canonical_data_selected_pval <- canonical_data_selected %>% select(c('pathway','group', 'log_p'))
canonical_data_selected_NES <- canonical_data_selected %>% select(c('pathway','group', 'NES'))


canonical_data_selected_pval_mt <- canonical_data_selected_pval %>% 
  as_tibble() %>% pivot_wider(id_cols = pathway,
                              names_from = group,
                              values_from = log_p,
                              values_fill = 0)%>%
  column_to_rownames('pathway')

canonical_data_selected_NES_mt <- canonical_data_selected_NES %>% 
  as_tibble() %>% pivot_wider(id_cols = pathway,
                              names_from = group,
                              values_from = NES,
                              values_fill = 0)%>%
  column_to_rownames('pathway')

#Change Column order
mt_colum_order <- c( "M0","Ribosome-Intermediate","HSP-Intermediate","Cytokine-Intermediate","Interferon","Cycling-G2M","Cycling-S","MGnD","MGnD-Antigen")
canonical_data_selected_NES_mt <- canonical_data_selected_NES_mt[,mt_colum_order]
canonical_data_selected_pval_mt <- canonical_data_selected_pval_mt[,mt_colum_order]

mt_row_order <- read.xlsx(paste0('results/canonical_data_total.xlsx'), sheet = 3)%>% pull(ordered_pathways)
canonical_data_selected_NES_mt <- canonical_data_selected_NES_mt[mt_row_order,]
canonical_data_selected_pval_mt <- canonical_data_selected_pval_mt[mt_row_order,]

#Edit rownames
rownames(canonical_data_selected_pval_mt) <- sub('REACTOME_','',rownames(canonical_data_selected_pval_mt))
rownames(canonical_data_selected_pval_mt) <- sub('BIOCARTA_','',rownames(canonical_data_selected_pval_mt))
rownames(canonical_data_selected_pval_mt) <- sub('WP_','',rownames(canonical_data_selected_pval_mt))
rownames(canonical_data_selected_pval_mt) <- gsub('_',' ',rownames(canonical_data_selected_pval_mt))

rownames(canonical_data_selected_NES_mt) <- sub('REACTOME_','',rownames(canonical_data_selected_NES_mt))
rownames(canonical_data_selected_NES_mt) <- sub('BIOCARTA_','',rownames(canonical_data_selected_NES_mt))
rownames(canonical_data_selected_NES_mt) <- sub('WP_','',rownames(canonical_data_selected_NES_mt))
rownames(canonical_data_selected_NES_mt) <- gsub('_',' ',rownames(canonical_data_selected_NES_mt))

rownames(canonical_data_selected_NES_mt) <- tolower(rownames(canonical_data_selected_NES_mt))
rownames(canonical_data_selected_pval_mt) <- tolower(rownames(canonical_data_selected_pval_mt))

canonical_data_selected_pval_mt <- as.matrix(canonical_data_selected_pval_mt)
canonical_data_selected_NES_mt <- as.matrix(canonical_data_selected_NES_mt)

#Set up legend
col_fun = colorRamp2(c(1.3, 2, 2.6), c("white", "orange", "red"))
col_fun(seq(1.3, 2.6))
lgd = Legend(col_fun = col_fun, title = "NES")

htmp_pval <- Heatmap(canonical_data_selected_NES_mt,
                     col = col_fun,
                     rect_gp = gpar(col = "black", lwd = 0.2),
                     cluster_columns = FALSE,
                     cluster_rows = FALSE,
                     width = ncol(canonical_data_selected_pval_mt)*unit(5, "mm"), 
                     height = nrow(canonical_data_selected_pval_mt)*unit(5, "mm"),
                     show_heatmap_legend = F,
                     cell_fun = function(j, i, x, y, w, h, fill) {
                       if(canonical_data_selected_pval_mt[i, j] > 2) {
                         grid.text("***", x, y)
                       } else if(canonical_data_selected_pval_mt[i, j] > 1.6) {
                         grid.text("**", x, y)
                       } else if(canonical_data_selected_pval_mt[i, j] > 1.3 ) {
                         grid.text("*", x, y)
                       }else if(canonical_data_selected_pval_mt[i, j] < 1.3) {
                         grid.text("", x, y)
                       }
                     }) 

pdf(file = paste0('plots/heatmap_canonical_pathways_', 'SubCellType', '.pdf'), width = 10, height = 12)
htmp_pval
draw(lgd, x = unit(0.5, "cm"), just = c("left"))
dev.off()


######## TRAJECTORY ASSESSMENT ####################################
Idents(dataset_filtered_reclust) <- 'Neutrophil_states_4'
DimPlot(dataset_filtered_reclust)
dataset_filtered_reclust_subset <- subset(dataset_filtered_reclust, idents = c('M0','MGnD','Cytokine-Intermediate',
                                                                             'Intermediate'))

dataset_filtered_reclust.slingshot <- RunSlingshot(srt = dataset_filtered_reclust, group.by = "Neutrophil_states_4", reduction = "umap", start = 'Immature')

# dataset_filtered_subset_slingshot@active.ident <- factor(dataset_filtered_subset_slingshot@active.ident,
#                                                         levels=c('M0','Ribosome-Intermediate','HSP-Intermediate',
#                                                                  'Interferon','MGnD-Antigen',
#                                                                  'MGnD','Cytokine-Intermediate'))
# levels(dataset_filtered_subset_slingshot)

Trajectory_plot <- ClassDimPlot(dataset_filtered_reclust.slingshot, group.by = "Neutrophil_states_4", reduction = "umap", lineages = paste0("Lineage", 1:2), lineages_span = 0.3,lineages_line_bg = "white", palette = "Paired")
Trajectory_plot


ggsave(
  paste0('cibersort/single_cell/Schulte_Schrepping/plots/','dataset_filtered_reclust.slingshot','Neutrophil_states_4' , '.pdf'),
  plot = Trajectory_plot,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 300,
  height = 300,
  units = c("mm"),
  dpi = 300)


#Dynamic Plot
dynamicplot <- DynamicPlot(
  srt = dataset_filtered_reclust.slingshot, lineages = c("Lineage1", "Lineage2"), group.by = "Neutrophil_states_4",
  exp_method = 'log1p',
  features = c('CD274','TGFB1', 'IL18R1','IL17RA'))

dynamicplot

ggsave(
  paste0('cibersort/single_cell/Schulte_Schrepping/plots/','dataset_filtered_reclust.slingshot','Neutrophil_states_4', '.pdf'),
  plot = dynamicplot,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 300,
  height = 300,
  units = c( "mm"),
  dpi = 300)

#Expression Plots
ExpDimPlot_trajectory <- ExpDimPlot(dataset_filtered_subset_xenon.slingshot, features = paste0("Lineage", 1:2), reduction = "UMAP", theme_use = "theme_blank")
ExpDimPlot_trajectory

ggsave(
  paste0('plots/ExpDimPlot_trajectory_','dataset_filtered_subset_slingshot_xenon' ,'.pdf'),
  plot = ExpDimPlot_trajectory,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 300,
  height = 300,
  units = c( "mm"),
  dpi = 300)







######## GRN analysis #######
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCENIC") 
packageVersion("SCENIC")
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather")




