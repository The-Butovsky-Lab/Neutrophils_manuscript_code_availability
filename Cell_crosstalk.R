#Load libraries
library(tidyverse)
library(dplyr)
library(openxlsx)
library(circlize)
library(graphics)
library(ComplexHeatmap)


setwd('/Users/kiliankleemann/Dropbox/Neutrophils Projects/Human_SK-5DC3')
#Database from NichenetR
#Cell Crosstalk 
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
write.xlsx(lr_network,paste0('lr_network.xlsx'))


#### NEUTRO MICROGLIA CROSSTALK ####
#Load datasets for analysis
crosstalk_prefix <- 'xxx'
data <- read.xlsx(paste0('results/',crosstalk_prefix, '.xlsx'))
data <- data %>% filter(log2FoldChange > 0.5 | log2FoldChange < -0.5)

sig_genes <- data %>% pull(gene)
sig_genes_network <- lr_network %>% filter(from %in% sig_genes)

#Import Interaction cell genes
interaction_celltype <- 'MICROGLIA_DEGS'

#Microglia
#DS Counts
geneExpression_Microglia <- read.xlsx("xxx")
geneExpression_Microglia_genes <- geneExpression_Microglia %>% filter(AVERAGE > 0.1) %>% pull(gene)

#Significant genes
DEGgeneExpression_Microglia <- read.xlsx("xxx")
DEGgeneExpression_Microglia_genes <- DEGgeneExpression_Microglia %>% filter(p_val< 0.001) %>% filter(pct.1 > 0.05) %>% pull(gene)
DEGgeneExpression_Microglia_data <- DEGgeneExpression_Microglia %>% filter(p_val< 0.001) %>% filter(pct.1 > 0.05)
DEGgeneExpression_Microglia_data$to <- DEGgeneExpression_Microglia_data$gene


#Final list 
sig_genes_network_final <- sig_genes_network %>% filter(to %in% DEGgeneExpression_Microglia_genes)
sig_genes_network_final$gene <- sig_genes_network_final$from
sig_genes_network_final <- merge(sig_genes_network_final,data,by=c("gene")) 
sig_genes_network_final$Receptor_Ligand <- paste0(sig_genes_network_final$from, "_",sig_genes_network_final$to)
sig_genes_network_final <- sig_genes_network_final[!duplicated(sig_genes_network_final$Receptor_Ligand),]
sig_genes_network_final <- merge(sig_genes_network_final,DEGgeneExpression_Microglia_data,by=c("to")) 


sig_genes_network_final <- sig_genes_network_final %>% group_by(gene.x) %>% mutate(Count = n())
sig_genes_network_final$Combined <- sig_genes_network_final$Count * sig_genes_network_final$log2FoldChange
sig_genes_network_final$FC_added <- sig_genes_network_final$log2FoldChange + sig_genes_network_final$avg_log2FC

#Final_filtering
sig_genes_network_final <- sig_genes_network_final %>% filter(FC_added > 0.5 | FC_added < -0.5)

#########---------Making Circlos Plots---------########
#Select Data
df <- as.data.frame(sig_genes_network_final)
df <- df   %>% 
  dplyr::select(c(`from`, `to`, `Combined`))
df <- df %>% arrange(desc(`Combined`))

#Setting Colour scheme
combined_max <- max(df$Combined)
combined_min <- min(df$Combined)
col_fun = colorRamp2(c(combined_min,0,combined_max), c("#112bcc","whitesmoke","orchid3"))
col_fun(seq(combined_min,combined_max, by = 0.01))

#Customize graphic parameters before initializing the circlos: 
circos.clear()
circos.par(canvas.xlim = c(-1.5, 1.5), 
           canvas.ylim = c(-1.5, 1.5),
           track.margin= c(0.01, 0.01),
           start.degree = -90,     #rotating circlos plot # of °
           #gap.degree = 0.8,
           "track.height" = 0.1)

#Assigning grid and annotation regions / size
chordDiagram(df, big.gap = 40, small.gap = 2,
             annotationTrack = c('grid','names'), 
             col=col_fun, 
             annotationTrackHeight = mm_h(2), 
             preAllocateTracks = list(track.height = mm_h(4)),h.ratio=0.4,
             transparency = 0.1,  
             directional = 1, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow",
             link.sort = TRUE, link.decreasing = TRUE,
             link.zindex = rank(df[[3]]))


#Assign Annotations 
circos.track(track.index = 1, panel.fun = function(x,y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1]+ mm_y(5), CELL_META$sector.index,
              cex = 2,
              facing = 'clockwise',
              niceFacing = T, 
              adj = c(0,0.9))
}, bg.border = NA)
title(str_wrap(paste0('Crosstalk_Neutrophils_', crosstalk_prefix, ' OLAH_MG_DEGs_expression_APOE34_vs_APOE33_AD_MCI_FC_Combined'),width =40))

#adding legend 
lgd_links = Legend(at = c(combined_min,0,combined_max), col_fun = col_fun, 
                   title_position = "topleft", title = "Regulatory Potential")
draw(lgd_links, x = unit(6, "mm"), y = unit(6, "mm"), just = c("left", "bottom"))

#highlighting sectors 
highlight.sector(df$from, track.index = 1, col = 'navyblue', text = 'NEUTROPHIL LIGANDS', cex = 0.7, text.col = 'white', facing = 'bending.inside', niceFacing = T)
highlight.sector(df$to, track.index = 1, col = 'cyan3', text = paste0(interaction_celltype,'RECEPTORS'), cex = 0.7, text.col = 'white', facing = 'bending.inside', niceFacing = T)

