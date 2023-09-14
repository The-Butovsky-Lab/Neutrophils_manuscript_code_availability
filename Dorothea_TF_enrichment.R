BiocManager::install('dorothea')
library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)
library(openxlsx)
library(ggrepel)

#Data and metadata
## We read the normalised counts and the experimental design 
Normalised_counts <- read.xlsx("results/xx/xxx")

#Experimental design
Experimental_design <- read.xlsx("xxx")


#Set experiment name
#file_postfix = "Spleen_F_12M_AD_APOE4_(remove_C09_merge_Het_and_HOM)_data"
file_postfix = "xxx"
#file_postfix = "Spleen_F_12M_AD_WT_data_APOE4_vs_APOE3_"
file_postfix
## We read the results from the differential analysis. 
DEGs <- read.xlsx(paste0("results/xxx/",file_postfix, "DEGene_statistics_pval05.xlsx")) %>% na.omit()

#Modification of layout
Normalised_counts_matrix <- Normalised_counts %>%
  remove_missing() %>%
  dplyr::mutate_if(~ any(is.na(.x)),~ if_else(is.na(.x),0,.x)) %>% 
  tibble::column_to_rownames(var = "gene") %>% 
  as.matrix() 

DEGs_matrix <- DEGs %>% 
  dplyr::select(gene, log2FoldChange) %>% 
  dplyr::filter(!is.na(log2FoldChange)) %>% 
  column_to_rownames(var = "gene")%>%
  as.matrix()


## We load Dorothea Regulons
data("dorothea_mm", package = "dorothea")
regulons <- dorothea_mm %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

###TFactivity TOP 25
tf_activities_stat <- dorothea::run_viper(DEGs_matrix, regulons,
                                          options =  list(minsize = 5, eset.filter = FALSE, 
                                                          cores = 1, verbose = FALSE, nes = TRUE))
tf_activities_stat_top20 <- tf_activities_stat %>%
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::rename(NES = "log2FoldChange") %>%
  dplyr::top_n(20, wt = abs(NES)) %>%
  dplyr::arrange(NES) %>% 
  dplyr::mutate(GeneID = factor(GeneID))

tf_activities_stat_top20_gg <- ggplot(tf_activities_stat_top20,aes(x = reorder(GeneID, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "cyan4", high = "orangered", 
                       mid = "white", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text( size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1, size =20, color = 'black'),
        axis.text.y = element_text(size =20,  color = 'black'),
        axis.line.y = element_line(color = 'black'),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Transcription Factors")


ggsave(paste0('plots/tf_activity/',file_postfix, 'cyan_orange.pdf'),
       plot = tf_activities_stat_top20_gg,
       device = NULL,
       path = NULL,
       #scale = 1,
       width = 20,
       height = 20,
       units = c("cm"),
       dpi = 600)


paletteLength <- 100
myColor <- colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

dorotheaBreaks <- c(seq(min(tf_activities_vector), 0, 
                        length.out=ceiling(paletteLength/2) + 1),
                    seq(max(tf_activities_vector)/paletteLength, 
                        max(tf_activities_vector), 
                        length.out=floor(paletteLength/2)))
dorothea_hmap <- pheatmap(tf_activities_counts_filter,
                          fontsize=14, 
                          fontsize_row = 8, 
                          fontsize_col = 8, 
                          color=myColor, 
                          breaks = dorotheaBreaks,
                          main = "Dorothea ABC", 
                          angle_col = 45,
                          treeheight_col = 0,  
                          border_color = NA,
                          cluster_rows = T)

