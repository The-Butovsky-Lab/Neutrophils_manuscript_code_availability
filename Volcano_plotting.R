#Loading packages
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')
BiocManager::install('DESeq2')
install.packages('tidyverse')
install.packages("openxlsx")

library(DESeq2)
library(EnhancedVolcano)
library(tidyverse)
library(openxlsx)
library(ggrepel)
library(grid)

######## SCATTER PLOT ####################
Volcano_directory <-  'xxx'
Volcano_directory_label <-  'xxx'
Volcano_prefix <- 'xxx'
volcano <- read.xlsx(paste0('results/',Volcano_directory,'statistics.xlsx'))
volcano$log_pval <- -log10(volcano$pvalue)
volcano$log_baseMean <- log10(volcano$baseMean)
sig_data <- volcano %>% filter(pvalue < 0.05)
nosig_data <- volcano %>% filter(pvalue > 0.05)

##Theme
theme_scatter <- theme(aspect.ratio = 1, 
                       panel.background = element_blank(),
                       panel.border=element_rect(fill=NA),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       strip.background=element_blank(),
                       axis.title  = element_text(size = 20, colour = "black"),
                       axis.text.x=element_text(colour="black",size =20),
                       axis.text.y=element_text(colour="black",size =20),
                       axis.ticks=element_line(colour="black"),
                       plot.margin=unit(c(1,1,1,1),"line"),
                       plot.title = element_text(colour="black",size =20))

###SELECT DATA - manual annotation -
select_UP <- read.xlsx(paste0('plots/volcanos/data_edit/',Volcano_directory_label,'DEGene_statistics_pval05.xlsx'),sheet = 2) %>% pull('UP')
select_DOWN <- read.xlsx(paste0('plots/volcanos/data_edit/',Volcano_directory_label,'DEGene_statistics_pval05.xlsx'),sheet = 2) %>% pull('DOWN')

#Select data
select_UP_data <- volcano %>% filter(gene %in% c(select_UP))
select_DOWN_data <- volcano %>% filter(gene %in% c(select_DOWN))

#Select dimensions
col_max <- max(sig_data$log2FoldChange)

scatter_volcano <- ggplot(volcano, aes(x= log2FoldChange, y= log_pval, label = gene)) +
  geom_point(data = sig_data , aes(fill = log2FoldChange),size = 3, shape = 21, color = 'black') +
  geom_point(data = nosig_data,size = 1, shape = 21, fill = 'grey')+
  scale_fill_gradient2(limits=c(-col_max, col_max), low="blue", mid="whitesmoke", high = "red", na.value = 'blue') +
  theme_scatter +
  # xlim(-1,1)+
  # ylim(0,9)+
  geom_vline(xintercept = 0.1,colour="grey", linetype = "longdash") +
  geom_vline(xintercept = -0.1,colour="grey", linetype = "longdash") +
  geom_hline(yintercept = 1.3, colour="grey", linetype = "longdash") +
  labs(title = str_wrap(paste0(Volcano_prefix),60), x = "Log2FC", y = "-log(pvalue)") +
  geom_text_repel(data = select_UP_data,
                  color = 'black', size = 6,fontface = 'italic',
                  min.segment.length = 0,
                   nudge_y = 3,
                  # nudge_x = 6 - select_UP_data$log2FoldChange,
                  nudge_x = 2,
                  max.overlaps = Inf,
                  box.padding = unit(0.5, 'mm'),
                  point.padding = unit(0.5, 'mm')) +
  geom_text_repel(data = select_DOWN_data,
                   color = 'black',size = 6, fontface = 'italic',
                   min.segment.length =0,
                    nudge_y = 3,
                   # nudge_x = -6 + abs(select_DOWN_data$log2FoldChange),
                   nudge_x = -2,
                   max.overlaps = Inf,
                   point.padding = unit(0.5, 'mm'))
scatter_volcano 


ggsave(paste0('plots/volcanos/',Volcano_prefix, '.pdf'),
       plot = scatter_volcano,
       device = NULL,
       path = NULL,
       #scale = 1,
       width = 20,
       height = 20,
       units = c("cm"),
       dpi = 600)

