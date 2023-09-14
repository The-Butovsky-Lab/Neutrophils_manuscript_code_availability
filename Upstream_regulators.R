#Comparison
#Load packages
list.of.packages <- c('tidyverse', 'DESeq2', 'openxlsx', 'RColorBrewer', 'VennDiagram',  'RColorBrewer',
                      'ggpubr','ggrepel', 'ggplot2', 'circlize', 'ComplexHeatmap','dplyr', 'plyr')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

#Load DEGs:
data1 <- read.xlsx(paste0('results/', 'xxx', '.xlsx'), sheet =4) 

#Upstream regulator gene sets (only pvalue < 0.05 and non 0 Z-score)
set1 <- na.omit(data1 %>% pull(xxx))
set2 <- na.omit(data1 %>% pull(xxx))
set3 <- na.omit(data1 %>% pull(xxx))

#### Identify overlap ####
IPA_upstream_overlap <- Reduce(intersect,list(set1,set2,set3))

#calculate overlap
phyper(307, 983, 13138-983, 3078,lower.tail=FALSE)


gene_set_comparison_name <- 'xxx'
gene_set_comparison <- 'xxx'


# Chart
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
venn <- venn.diagram(
  x = list(set1,set2,set3),
  category.names = c("xxx" , "xxx","xxx"),
  filename = paste0('plots/overlap_comparison/',gene_set_comparison_name,'.png'),
  output=TRUE,
  # Output features
  imagetype = "png",
  height = 900, 
  width = 900, 
  resolution = 600,
  compression = "lzw",
  #Circles
  lwd = 1,
  col=c("#440154ff", '#21908dff','orangered'), #,'blue'),
  fill = c(alpha("#440154ff",0.3),alpha('#21908dff',0.3),alpha('orangered',0.3)),#, alpha('blue',0.3)),
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  # Set names
  cat.cex = 0.5,
  cat.default.pos = "outer",
  cat.pos = c(8,-10,5),
  cat.dist = c(0.2,0.3,0.25),
  cat.fontface = "bold",
  cat.col = c("#440154ff",'#21908dff','orangered'),
  cat.fontfamily = "sans")


#########---------SET UP DATA ---------########
DEGs_1 <- read.xlsx(paste0('results/IPA/xxx'))
DEGs_2 <- read.xlsx(paste0('results/IPA/xxx'))
DEGs_3 <- read.xlsx(paste0('results/IPA/xxx'))

DEGs_1$Group <- 'xxx'
DEGs_2$Group <- 'xxx'
DEGs_3$Group <- 'xxx'


barplot_data <- rbind(DEGs_1,DEGs_2,DEGs_3)
barplot_data

selected_regulators <-  c('xxx')
  
barplot_data_overlap <- barplot_data %>% filter(Upstream.Regulator %in% selected_regulators)
barplot_data_overlap <- barplot_data_overlap %>% filter(`p-value.of.overlap` < 0.05)
barplot_data_overlap$zscore <- as.numeric(barplot_data_overlap$`Activation.z-score`)*1
write.xlsx(barplot_data_overlap, paste0('results/xxx.xlsx'))

#########---------Making Barplots ---------########
### SET THEME ###
lineWidth = 1
pointSize = 30
theme_barplot <- theme(
  text = element_text(size = pointSize, colour = "black"),
  rect = element_blank(),
  line = element_line(size = lineWidth, colour = "black"),
  plot.title  = element_text(color="black", size=1, face="bold.italic"),
  axis.title  = element_text(size = pointSize, colour = "black"),
  axis.text.x  = element_text(size = pointSize , colour = "black", angle = 45, hjust = 1, face = 'italic'),
  axis.text.y  = element_text(size = pointSize , colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # legend.title = element_blank(),
  legend.position="right",
  legend.text = element_text(size = pointSize, colour = "black"),
  legend.key.height = unit(0.2, "cm"),
  legend.key.width = unit(0.2, "cm"),
  axis.line = element_line(size = lineWidth, colour = "black"),
  axis.ticks = element_line(size = lineWidth, colour = "black"),
  plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

######
barplot <- ggplot(barplot_data_overlap, aes(x = reorder(Upstream.Regulator, zscore), y = zscore, fill = Group)) + 
  geom_bar(position = 'dodge', stat = 'summary',width = 0.7, colour="black", alpha = 0.8) +
  expand_limits(x = 0, y = 0) +
  theme(panel.spacing = unit(1, "lines")) +
  labs(x = NULL, y = c("Log2FC")) +
  geom_hline(yintercept = 0, colour="black") +
  ylim(-3,3) +
  theme_barplot +
  scale_y_continuous(expand = c(0, 0.1, .2, 0))
barplot

dev.off()
pdf(file = paste0('plots/overlap_comparison/Barplot_overlap_',gene_set_comparison_name,'.pdf'), pointsize = 10, width = 30,height =10)
barplot
dev.off()


#### Trajectory plot Upstream_regulators ####
barplot_data_overlap_filter <- read.xlsx(paste0('results/IPA/xxx.xlsx'))
barplot_data_overlap_filter <- barplot_data_overlap_filter %>% select(c('Upstream.Regulator','Group','zscore','p-value.of.overlap','Molecule.Type'))


selected_regulators_group1 <-  c('xxx')
selected_regulators_group2 <- c('xxx')

barplot_data_overlap_filter_extra <- barplot_data_overlap_filter %>%filter(Upstream.Regulator %in% selected_regulators_group2)

#Reorderlevels
barplot_data_overlap_filter_extra$Group_1 <- factor(barplot_data_overlap_filter_extra$Group_1, levels = c('xxx'))

#trajectoryPlot
lineWidth = 1
pointSize = 20
Upstream_regulator_trajectory <- ggplot(barplot_data_overlap_filter, aes(x=Group_1, y = zscore, group = Upstream.Regulator), fill = Upstream.Regulator)+
  geom_hline(yintercept = 0, color = 'grey')+
  #geom_line(aes(color = PROTEIN), linetype="dashed", size=1) +
  geom_smooth(method = loess,aes(color = Upstream.Regulator), size=2) +
  geom_point(aes(color = Upstream.Regulator),size = 6) +
  theme(text = element_text(size = pointSize, colour = "black"),
        rect = element_blank(),
        line = element_line(size = lineWidth, colour = "black"),
        plot.title  = element_text(color="black", size = pointSize),
        plot.background=element_rect(fill="transparent",colour=NA),
        axis.title.x = element_text(color="black", size = pointSize),
        axis.title.y = element_text(color="black", size = pointSize),
        axis.text.x  = element_text(size = pointSize , colour = "black", angle = 90, vjust=0.1),
        axis.text.y  = element_text(size = pointSize , colour = "black"),
        axis.ticks = element_line(size = lineWidth, colour = "black"),
        axis.ticks.length=unit(.25, "cm"),
        axis.line = element_line(size = lineWidth, colour = "black"),
        legend.position = "right",
        legend.title = element_text(size = pointSize , colour = "black"),
        legend.text = element_text(size = pointSize , colour = "black"),
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.background = element_blank(),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        panel.background = element_rect(fill="transparent",colour=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(size = lineWidth, colour = "grey", linetype = 'dashed'))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  ggtitle(paste0('Zscore of upstream regulators ',' across Conditions'))+
  scale_y_continuous(expand = c(0, 0.1, .05, 0)) 
  # scale_color_manual(values = c(cols)) +
  # scale_fill_manual(values = c(cols)) 
Upstream_regulator_trajectory


ggsave(paste0('plots/trajectory/','xxx','.pdf'),
       Upstream_regulator_trajectory,
       dpi = 300)

write.xlsx(barplot_data_overlap_filter_extra_heatmap, paste0('results/IPA/xxx.xlsx'))

#### Heatmap Upstream Regulators ####
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

barplot_data_overlap_filter_extra_heatmap <- pivot_wider(barplot_data_overlap_filter,
                                                         names_from = Group_1,
                                                         values_from = zscore,
                                                         id_cols = Upstream.Regulator)

barplot_data_overlap_filter_extra_heatmap <- barplot_data_overlap_filter_extra_heatmap %>% column_to_rownames('Upstream.Regulator')

reordered_index <- c(
  grep("xxx", names(barplot_data_overlap_filter_extra_heatmap), ignore.case = T),
  grep("xxx", names(barplot_data_overlap_filter_extra_heatmap), ignore.case = T),
  grep("xxx", names(barplot_data_overlap_filter_extra_heatmap), ignore.case = T),
  grep("xxx", names(barplot_data_overlap_filter_extra_heatmap), ignore.case = T))

reordered_data <- barplot_data_overlap_filter_extra_heatmap %>%
  select(reordered_index)

data_norm <- t(apply(reordered_data, 1, cal_z_score))

pdf(file = paste0('plots/heatmaps/', 'xxx', '.pdf'), pointsize = 10)
phm_full <- pheatmap(data_norm,
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = 2,
                     cluster_cols = F,
                     #cutree_cols = 4,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 10,
                     #cellwidth = 25,
                     scale = 'row'
                     #cellheight = 3,
                     )
phm_full
dev.off()


#### Cumulative Upstream regulators ####
#Adding Zscores for comparisons together
Cumulative_upstream_regulator <- read.xlsx('results/IPA/xxx')

select_molecules <- read.xlsx('results/IPA/xxx', sheet =2)



regulators_down<-  c('xxx')
regulators_up <- c('xxx')

#Selecting classes of regulators 
regulator_class <- 'transporter'
selected_regulators <- select_molecules %>% filter(Molecule_type == regulator_class) %>% pull(Upstream.Regulator)

selected_regulators <- regulators_up
#Set colors
cols <- brewer.pal(9, "BuPu")
cols <- cols[4:5]

color_chosen <- c('darkorchid1','darkorchid4')

#Select subset
Cumulative_upstream_regulator_selected <- Cumulative_upstream_regulator %>% filter(Upstream.Regulator %in% selected_regulators)
Cumulative_upstream_regulator_longer <- Cumulative_upstream_regulator_selected %>% pivot_longer(cols = -(1:1),
                                                                                values_to = c("Cumulative_zscore"),
                                                                                names_to = c("Condition"))
Cumulative_upstream_regulator_longer$Condition <- factor(Cumulative_upstream_regulator_longer$Condition, levels = c('xxx'))

#Make plot
Upstream_regulator_trajectory <- ggplot(Cumulative_upstream_regulator_longer, aes(x=Condition, y = Cumulative_zscore, group = Upstream.Regulator, fill = Upstream.Regulator))+
  geom_hline(yintercept = 0, color = 'grey')+
  geom_smooth(method = loess,aes(color = Upstream.Regulator), size=2) +
  geom_point(aes(fill = Upstream.Regulator),color='black',size = 6, pch=21) +
  theme(text = element_text(size = pointSize, colour = "black"),
        rect = element_blank(),
        line = element_line(size = lineWidth, colour = "black"),
        plot.title  = element_text(color="black", size = pointSize),
        plot.background=element_rect(fill="transparent",colour=NA),
        axis.title.x = element_text(color="black", size = pointSize),
        axis.title.y = element_text(color="black", size = pointSize),
        axis.text.x  = element_text(size = pointSize , colour = "black", angle = 90, vjust=0.1),
        axis.text.y  = element_text(size = pointSize , colour = "black"),
        axis.ticks = element_line(size = lineWidth, colour = "black"),
        axis.ticks.length=unit(.25, "cm"),
        axis.line = element_line(size = lineWidth, colour = "black"),
        legend.position = "right",
        legend.title = element_text(size = pointSize , colour = "black"),
        legend.text = element_text(size = pointSize , colour = "black"),
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.background = element_blank(),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        panel.background = element_rect(fill="transparent",colour=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(size = lineWidth, colour = "grey", linetype = 'dashed'))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  ggtitle(str_wrap(paste0('Zscore of ',regulator_class,' upstream regulators ',' across Conditions'),40))+
  scale_y_continuous(expand = c(0, 0.5, .1, 0)) +
  scale_color_manual(values = c(color_chosen)) +
  scale_fill_manual(values = c(color_chosen))
Upstream_regulator_trajectory


ggsave(paste0('plots/trajectory/','xxx', regulator_class,'.pdf'),
       Upstream_regulator_trajectory,
       dpi = 300)



#### Dotplot for Upstream regulators ####
barplot_data_overlap_filter <- read.xlsx(paste0('results/IPA/xxx'))
barplot_data_overlap_filter <- barplot_data_overlap_filter %>% select(c('Upstream.Regulator','Group','zscore','p-value.of.overlap','Molecule.Type'))

barplot_data_overlap_filter$log_pval <- -log10(barplot_data_overlap_filter$`p-value.of.overlap`)
barplot_data_overlap_filter$Group <- factor(barplot_data_overlap_filter$Group, levels = c('xxx','xxx','xxx'))

dotplot <- ggplot(barplot_data_overlap_filter,aes(x = Group, y = reorder(Upstream.Regulator,zscore))) + 
  geom_point(aes(size = log_pval),color ='black', stroke = 4) +
  geom_point(aes(size = log_pval, color = zscore), stroke = 3) +
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
        #panel.margin=unit(.5, "lines"),
        #panel.border = element_rect(color = "black", fill = NA, size = 1),
        #panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0,size = pointSize),
        #strip.background = element_rect(size = 1))+
        strip.text.y = element_text(size=10,angle = 90)) +
  #facet_wrap(~group_pathway, strip.position = "left", scales = "free_y") +
  scale_color_gradient2(limits=c(-4,4), low="cyan4", mid="white", high = "orangered", na.value = 'orangered') +
  #scale_fill_gradientn(colours = c("darkred","red2", "red", "orange","white","darkblue"),
  #                       values = c(1,0.8,0.6,0.4,0.2,0)) +
  #scale_colour_manual(values = c('navyblue','antiquewhite','firebrick3'), limits = c("-1", "0", "3"))+
  labs(title = str_wrap(paste0("Upstream regulators"), 40))+
  #facet_grid(~Molecule.Type, scales = "free_x", space='free', switch = "both")+ # for x-axis
  facet_grid(Molecule.Type~., scales = "free_y", space='free', switch = "both")
  #coord_flip()
#facet_grid(.~group_pathway, scales = "free", switch = "y", space = "free_x") + 
#geom_text(position = position_dodge(width = -5), aes(y=group_pathway, x=0))
#theme(strip.placement = "outside")
dotplot

dev.off()
pdf(file = paste0('plots/pathways/IPA/Dotplot_','xxx', '.pdf'), width = 10, height = 20)
dotplot
dev.off()

ggsave(paste0('plots/pathways/IPA/Dotplot_','xxx', '.pdf'),
       dotplot,
       dpi = 600)




#Scatter plot positive and negative genes in overlap 
scatter_data <- merge(DEG_1,DEG_2, by = 'gene')

theme_scatter <- theme(aspect.ratio = 1, 
                       panel.background = element_blank(),
                       panel.border=element_rect(fill=NA),
                       panel.grid.major = element_line(linetype = "dashed"),
                       panel.grid.minor = element_blank(),
                       strip.background=element_blank(),
                       axis.text.x=element_text(colour="black"),
                       axis.text.y=element_text(colour="black"),
                       axis.ticks=element_line(colour="black"),
                       plot.margin=unit(c(1,1,1,1),"line"))

scatter <- ggplot(scatter_data, aes(x= log2FoldChange.x, y= log2FoldChange.y, label = gene)) +
  geom_smooth(data = scatter_data, method = "lm", se=TRUE, color="red", formula = y ~ x) +
  stat_cor(data = scatter_data, color = 'red', label.x.npc = "left", label.y.npc = "top") +
  #geom_smooth(data = M0_data, method = "lm", se=TRUE, color="dodgerblue", formula = y ~ x) +
  #stat_cor(data = M0_data, color = 'dodgerblue', label.x.npc = "center", label.y.npc = "bottom") +
  geom_point(size=0.5, color = 'black') +
  #geom_point(data = DEG_data, color = 'black', size = 0.5)+
  # geom_point(data = IFNg_data, color = 'purple', size = 1) +
  # geom_point(data = Antigen_data, color = 'chartreuse4', size = 1)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  theme +
  labs(title = paste0('Comparison','_',name_13,"_",name_15),
       x = paste0("Log2FC_",name_13),
       y = paste0("Log2FC_",name_15)) 
#geom_label(data = IFNg_data, color = 'orange', size = 2, label.padding = unit(0.1,"lines")) +
#geom_label(data = Antigen_data, color = 'chartreuse4', size = 2, label.padding = unit(0.1, "lines")) 
# geom_label_repel(data = IFNg_data, 
#                  color = 'purple', size = 3, 
#                  min.segment.length = 0,
#                  max.overlaps = Inf,
#                  label.padding = unit(0.1,"lines")) +
# geom_label_repel(data = Antigen_data, 
#                  color = 'chartreuse4',size = 3, 
#                  min.segment.length = 0,
#                  max.overlaps = Inf,
#                  label.padding = unit(0.1,"lines"))
scatter 


dev.off()
pdf(file = paste('plots/overlap_comparison/Scatterplot_overlap_',gene_set_comparison_name,'.pdf'), pointsize = 10, width =6, height =6)
scatter
dev.off()





