# References: https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html

## Gene-level differential expression analysis using DESeq2
## !requireNamespace is used to check whether the package is installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dendextend", force = TRUE)
BiocManager::install("fgsea", force = TRUE)
BiocManager::install("tidyverse", force = TRUE)
BiocManager::install("enrichplot", force = TRUE)
BiocManager::install("org.Mm.eg.db", force = TRUE)
BiocManager::install("org.Hs.eg.db", force = TRUE)

library(BiocGenerics)
library(tximport)
library(S4Vectors)
library(DESeq2)
library(openxlsx)
library(biomaRt)
library(tidyverse)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(cowplot)
library(RColorBrewer)
library(sva)
library(ggpubr)
library(ggrepel)

#Set Working directory 
setwd("~/Dropbox/Neutrophils Projects/Human_SK-5DC3/")
meta <- readxl::read_excel(paste0(workdir,expdir,"metadata.xlsx"))
quantfiles <- file.path(workdir,expdir,"transcripts_quant/", meta$Well, "quant.sf")

##---------------------------------FILE PREPARATION----------------------------------##
# ## useMart is a function used to connect to the selected BioMart database and dataset (ListEnsemblArchives() to list host sites)
# ensembl_hs_mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# ## getBM is a function used to retrieve information from the BioMart database
# ensembl_df <- getBM(attributes = c("ensembl_transcript_id_version", "ensembl_gene_id_version", "external_gene_name"), mart = ensembl_hs_mart) %>%
#   as_tibble() %>%
#   mutate(noVersion = sub('\\..*', '', ensembl_transcript_id_version))%>%
#   mutate(noVersion_ENSG = sub('\\..*', '', ensembl_gene_id_version))
# 
# # Build tx2gene table from ensembl_df
# tx2gene <- ensembl_df %>%
#   select(noVersion, external_gene_name) %>%
#   rename(c('Transcript' = "noVersion" , "Gene" = "external_gene_name" )) %>%
#   data.frame()
# 
# tx2gene <- read.xlsx('tx2gene.xlsx')

#Create a metadata 
metadata_MCSA <- read.xlsx("Public_data/MCSA_blood/metadata/MCSA_individual_human_metadata.xlsx", sheet = 1)
immune_dev_res_samples_neutro_score_greater70_samples <- read.xlsx("Public_data/MCSA_blood/immune_dev_res_samples_neutrophil_score_greater70_samples.xlsx", sheet = 1) %>% pull('immune_dev_res_samples_neutro_score_greater70_samples')

#Samples with more than 70% predited neutrophil content
CN_MCI_data <- metadata_MCSA %>% filter(diagnosis %in% c('CN','MCI'))
CN_MCI_neutrophil_70_data <- CN_MCI_data %>% filter(individualID %in% immune_dev_res_samples_neutro_score_greater70_samples)

CN_MCI_data_samples <- CN_MCI_data %>% pull(ID_specific_extra)
CN_MCI_neutrophil_70_data_samples <- CN_MCI_neutrophil_70_data %>% pull(ID_specific_extra)

#Subset metadata with only 70% neutrophil data
#Female
Female_data <- CN_MCI_neutrophil_70_data %>% filter(sex == "female")
Female_APOE33_data <- CN_MCI_neutrophil_70_data %>% filter(apoeGenotype == c('APOE_33') & sex == c('female')) 
Female_APOE34_data <- CN_MCI_neutrophil_70_data %>% filter(apoeGenotype == c('APOE_34') & sex == c('female')) 
Female_CN_data_APOE34_vs_APOE33 <- CN_MCI_neutrophil_70_data %>% filter(diagnosis == c('CN') & sex == c('female')& apoeGenotype == c('APOE_34','APOE_33')) 
Female_CN_data_APOE34_vs_APOE23 <- CN_MCI_neutrophil_70_data %>% filter(diagnosis == c('CN') & sex == c('female')& apoeGenotype == c('APOE_34','APOE_23')) 

Female_MCI_data_APOE34_vs_APOE33 <- CN_MCI_neutrophil_70_data %>% filter(diagnosis == c('MCI') & sex == c('female')& apoeGenotype == c('APOE_34','APOE_33')) 
Female_MCI_data_APOE34_vs_APOE23 <- CN_MCI_neutrophil_70_data %>% filter(diagnosis == c('MCI') & sex == c('female')& apoeGenotype == c('APOE_34','APOE_23')) 

#Male
Male_APOE33_data <- CN_MCI_neutrophil_70_data %>% filter(apoeGenotype == c('APOE_33') & sex == c('male')) 
Male_APOE34_data <- CN_MCI_neutrophil_70_data %>% filter(apoeGenotype == c('APOE_34') & sex == c('male')) 

Male_CN_data_APOE34_vs_APOE33 <- CN_MCI_neutrophil_70_data %>% filter(diagnosis == c('CN') & sex == c('male')& apoeGenotype == c('APOE_34','APOE_33')) 
Male_MCI_data_APOE34_vs_APOE33 <- CN_MCI_neutrophil_70_data %>% filter(diagnosis == c('MCI') & sex == c('male')& apoeGenotype == c('APOE_34','APOE_33')) 

#EXTRACT SAMPLE NAMES
#Female
Female_data_samples  <- Female_data %>% pull(ID_specific_extra)
Female_APOE33_data_samples <- Female_APOE33_data %>% pull(ID_specific_extra)
Female_APOE34_data_samples <- Female_APOE34_data %>% pull(ID_specific_extra)
Female_CN_data_APOE34_vs_APOE33_samples <- Female_CN_data_APOE34_vs_APOE33 %>% pull(ID_specific_extra)

Female_CN_data_APOE34_vs_APOE23_samples <- Female_CN_data_APOE34_vs_APOE23 %>% pull(ID_specific_extra)

Female_MCI_data_APOE34_vs_APOE33_samples <- Female_MCI_data_APOE34_vs_APOE33 %>% pull(ID_specific_extra)
Female_MCI_data_APOE34_vs_APOE23_samples <- Female_MCI_data_APOE34_vs_APOE23 %>% pull(ID_specific_extra)

#Male
Male_APOE33_data_samples <- Male_APOE33_data %>% pull(ID_specific_extra)
Male_APOE34_data_samples <- Male_APOE34_data %>% pull(ID_specific_extra)
Male_CN_data_APOE34_vs_APOE33_samples <- Male_CN_data_APOE34_vs_APOE33 %>% pull(ID_specific_extra)
Male_MCI_data_APOE34_vs_APOE33_samples <- Male_MCI_data_APOE34_vs_APOE33 %>% pull(ID_specific_extra)


#Expression data 
expression_matrix_MCSA <- read.xlsx('Public_data/MCSA_blood/data/mcsa_rnaseq_rawcount_edit_names.xlsx')
# expression_matrix_CN_MCI_data <- select(expression_matrix_MCSA,matches(c("GeneName", CN_MCI_data_samples)))
# expression_matrix_CN_MCI_data <- expression_matrix_CN_MCI_data %>% distinct(GeneName,.keep_all = TRUE)
# expression_matrix_CN_MCI_data <- expression_matrix_CN_MCI_data %>% column_to_rownames('GeneName')

expression_matrix_CN_MCI_data_neutrophil_70 <- select(expression_matrix_MCSA,matches(c("GeneName",CN_MCI_neutrophil_70_data_samples)))

#Subset expression data for females
expression_matrix_Female_data <- select(expression_matrix_CN_MCI_data_neutrophil_70,c("GeneName",Female_data_samples))

expression_matrix_Female_APOE33_data <- dplyr::select(expression_matrix_CN_MCI_data_neutrophil_70,matches(c("GeneName",Female_APOE33_data_samples)))
expression_matrix_Female_APOE34_data <- select(expression_matrix_CN_MCI_data_neutrophil_70,matches(c("GeneName",Female_APOE34_data_samples)))

expression_matrix_Female_CN_data_APOE34_vs_APOE33 <- select(expression_matrix_CN_MCI_data_neutrophil_70,matches(c("GeneName",Female_CN_data_APOE34_vs_APOE33_samples))) 
expression_matrix_Female_MCI_data_APOE34_vs_APOE33 <- select(expression_matrix_CN_MCI_data_neutrophil_70,matches(c("GeneName",Female_MCI_data_APOE34_vs_APOE33_samples))) 

expression_matrix_Female_CN_data_APOE34_vs_APOE23 <- select(expression_matrix_CN_MCI_data_neutrophil_70,matches(c("GeneName",Female_CN_data_APOE34_vs_APOE23_samples))) 

#Subset expression data for males
expression_matrix_Male_APOE33_data <- select(expression_matrix_MCSA,matches(c("GeneName",Male_APOE33_data_samples)))
expression_matrix_Male_CN_data_APOE34_vs_APOE33 <- select(expression_matrix_CN_MCI_data,matches(c("GeneName",Male_CN_data_APOE34_vs_APOE33_samples))) 
expression_matrix_Male_MCI_data_APOE34_vs_APOE33 <- select(expression_matrix_CN_MCI_data,matches(c("GeneName",Male_MCI_data_APOE34_vs_APOE33_samples))) 


#Analyzis parameters
file_prefix = 'Female_APOE33_data_CHECK_'
data = expression_matrix_Female_APOE33_data
experiment = Female_APOE33_data
experiment


#### Explore filtered metadata ####
# lineWidth = 1
# pointSize = 10
# histogram <- ggplot(experiment, aes(x=Age)) + 
#   geom_histogram(aes(y=..density..), color="black", fill="lightblue", linetype="dashed", binwidth = 4)+
#   geom_density(alpha=.2, fill="cyan",adjust = 2) + 
#   theme(text = element_text(size = pointSize, colour = "black"),
#         rect = element_blank(),
#         line = element_line(size = lineWidth, colour = "black"),
#         plot.title  = element_text(color="black", size=20),
#         axis.title.x = element_text(size = pointSize , colour = "black"),
#         axis.title.y = element_text(size = pointSize , colour = "black"),
#         axis.text.x  = element_text(size = pointSize , colour = "black"),
#         axis.text.y  = element_text(size = pointSize , colour = "black"),
#         axis.ticks = element_line(size = lineWidth, colour = "black"),
#         axis.line = element_line(size = lineWidth, colour = "black"),
#         legend.position = "right",
#         legend.title = element_text(size = pointSize , colour = "black"),
#         legend.text = element_text(size = pointSize , colour = "black"),
#         legend.key.height = unit(1, "cm"),
#         legend.key.width = unit(0.5, "cm"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.margin=unit(.5, "lines"))+
#   scale_y_continuous(expand = c(0, 0, .05, 0)) +
#   scale_x_continuous(expand = c(0, 0, 0, 0), limits=c(0,90)) 
# histogram
# 
# 
# pdf(file = paste0('plots/histogram/','Histogram_age','.pdf'), pointsize = 10, width = 25, height= 5)
# print(histogram)
# dev.off()


####
#Filtering GeneNames
data <- data %>% rownames_to_column('GeneName')
data_filtered <- data %>% distinct(GeneName,.keep_all = TRUE)

#Filter minimum counts
data_filtered <- data_filtered %>% column_to_rownames('GeneName')
idx <- rowSums(data_filtered) >= 20
data_filtered <- data_filtered[idx,]






##### DESEQ2 NALYSIS BEGINS HERE #####
## Create DESeq2Dataset object with variables of interest.
dds <- DESeqDataSetFromMatrix(data_filtered, colData = experiment, design = ~  diagnosis)  

#prefiltering on minimum of 5 reads (NOT REQUIRED AS DESEQ2 OPTIMIZES AND FILTERS AUTOMATICALLY)
dds <- estimateSizeFactors(dds)
# idx <- rowSums(counts(dds, normalized = TRUE)) >= 30
# dds <- dds[idx,]

#Setting the reference level (control group to compare against)
dds@colData@listData$diagnosis <- relevel(dds@colData@listData$diagnosis, ref = "CN")
#dds@colData@listData$clinical_2 <- relevel(dds@colData@listData$clinical_2, ref = "MCI")
#dds@colData@listData$apoeGenotype <- relevel(dds@colData@listData$apoeGenotype, ref = "APOE_33")
#dds@colData@listData$sex <- relevel(dds@colData@listData$sex, ref = "F")


## Run DESeq analysis to gather differential expression results
#Run DESeq (LRT)
dds_run <- DESeq(dds, test = 'LRT', reduced = ~ 1)
# Run DESeq (Wald)
#dds_run <- DESeq(dds, betaPrior = F)

# View names of estimated effects
resultsNames(dds_run)

# Create table of effects showing log2fold, standard error, test stats, and p-vals
# LRT
dds_result <- results(dds_run, pAdjustMethod = 'BH', independentFiltering = T, cooksCutoff = T)

# Wald
# dds_result <- results(dds_run, contrast = c('diagnosis','MCI', 'CN'), independentFiltering = T, pAdjustMethod = 'BH')
# dds_result <- lfcShrink(dds_run, contrast = c('diagnosis', 'MCI', 'CN'), res = dds_result, type = 'normal')

dds_result
summary(dds_result, alpha = 0.05)


## View Total number of normalized counts per sample
dds_run <- estimateSizeFactors(dds_run)
DS_norm_counts <- counts(dds_run, normalized = TRUE)
colSums(DS_norm_counts)

## Plot dispersion estimates
#plotDispEsts(dds_run)

#plotCounts(dds_run, 'Fus', intgroup = 'condition')

# Transform counts for data visualization
vst <- vst(dds_run, blind=TRUE)

# Plot PCA 
# Add nametags
z <- plotPCA(vst, intgroup=c('diagnosis'), ntop = 300)
#z + geom_label(aes(label = experiment$position))
theme<-theme(aspect.ratio = 1, panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
pdf(file = paste0('Public_data/MCSA_blood/plots/PCA/', file_prefix,'', '.pdf'), pointsize = 10)
z + geom_point() +
  #geom_label(aes(label = CN_MCI_data$individualID)) +
  theme + 
  scale_color_brewer(palette = 'Set1') 
dev.off()


pdf(file = paste0('plots/PCA/', file_prefix,'Ellipse', '.pdf'), pointsize = 10)
z + geom_point() + 
  stat_ellipse(aes(fill=group), geom = "polygon", alpha=0.5,level = 0.5) +
  theme +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.text = element_text(size=20)) +
  scale_color_brewer(palette = 'Paired') +
  scale_fill_brewer(palette = 'Paired')
dev.off()

#Filter unwated genes
#unwanted_genes = paste(c('^Gm', '^mt-', '^Vmn', '^Rpl', '^Rps', '^Olfr','Rik$'), collapse = '|')

# Build significant gene table and extract list of sorted DE genes.
sig_res <- dds_result %>%
  data.frame() %>%
  rownames_to_column(var='gene_name') %>%    # makes a 'gene' column using column1
  as_tibble %>%
  filter(pvalue < 0.05) %>%
  arrange(pvalue) #%>% 
#filter(!str_detect(gene, unwanted_genes))

sorted_DEGenes <- sig_res$gene_name

# Export significant gene count data
name_list <- c('gene_name', (paste0(dds_run$individualID, "_", dds_run$diagnosis, "_",dds_run$apoeGenotype, "_", dds_run$sex)))
name_list                        

sig_export <- DS_norm_counts %>%
  data.frame() %>%
  rownames_to_column(var = 'gene_name') %>%
  as_tibble() %>%
  `colnames<-`(name_list) %>%
  semi_join(sig_res) %>%
  dplyr::slice(match(sorted_DEGenes, gene_name)) #%>%
#filter(!str_detect(gene, unwanted_genes)) 

# Column Names
check_cols <- c(colnames(sig_export))
check_cols

# #Manual Reordering of Columns (if necessary)
reordered_index <- c(
  grep("CN_APOE_22_m", names(sig_export), ignore.case = T),
  grep("CN_APOE_23_m", names(sig_export), ignore.case = T),
  grep("CN_APOE_33_m", names(sig_export), ignore.case = T),
  grep("CN_APOE_24_m", names(sig_export), ignore.case = T),
  grep("CN_APOE_34_m", names(sig_export), ignore.case = T),
  grep("CN_APOE_44_m", names(sig_export), ignore.case = T),
  
  grep("CN_APOE_22_f", names(sig_export), ignore.case = T),
  grep("CN_APOE_23_f", names(sig_export), ignore.case = T),
  grep("CN_APOE_33_f", names(sig_export), ignore.case = T),
  grep("CN_APOE_24_f", names(sig_export), ignore.case = T),
  grep("CN_APOE_34_f", names(sig_export), ignore.case = T),
  grep("CN_APOE_44_f", names(sig_export), ignore.case = T),
  
  grep("MCI_APOE_22_m", names(sig_export), ignore.case = T),
  grep("MCI_APOE_23_m", names(sig_export), ignore.case = T),
  grep("MCI_APOE_33_m", names(sig_export), ignore.case = T),
  grep("MCI_APOE_24_m", names(sig_export), ignore.case = T),
  grep("MCI_APOE_34_m", names(sig_export), ignore.case = T),
  grep("MCI_APOE_44_m", names(sig_export), ignore.case = T),
  
  grep("MCI_APOE_22_f", names(sig_export), ignore.case = T),
  grep("MCI_APOE_23_f", names(sig_export), ignore.case = T),
  grep("MCI_APOE_33_f", names(sig_export), ignore.case = T),
  grep("MCI_APOE_24_f", names(sig_export), ignore.case = T),
  grep("MCI_APOE_34_f", names(sig_export), ignore.case = T),
  grep("MCI_APOE_44_f", names(sig_export), ignore.case = T)
)

# 
sig_export <- sig_export %>%
  dplyr::select(c(1, reordered_index))

check_cols <- c(colnames(sig_export))
check_cols

# Normalized counts
check_counts <- DS_norm_counts %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  dplyr::select(c(1, reordered_index)) %>%
  as_tibble %>% 
  `colnames<-`(check_cols) #%>%
#filter(!str_detect(gene, unwanted_genes))


# Results
check_res <- dds_result %>%
  as.data.frame() %>%
  rownames_to_column('gene') #%>%
#filter(!str_detect(gene, unwanted_genes))

#TPM File
# TPM_name_list <- c((paste0(dds_run$position, "_", dds_run$clinical_2, "_",dds_run$apoe, "_", dds_run$sex)))
# TPM_name_list  
# TPM  <- txicounts %>%
#   as.data.frame() %>%
#   `colnames<-`(TPM_name_list) 
# 

# Output files
#Statistics
write.xlsx(sig_res, file = paste0('Public_data/MCSA_blood/results/', file_prefix, 'DEGene_statistics_pval05.xlsx'), overwrite = T)
#write.xlsx(sig_res_padj, file = paste0('results/', file_prefix, 'DEGene_statistics_padj05.xlsx'), overwrite = T)
write.xlsx(check_res, file = paste0('Public_data/MCSA_blood/results/', file_prefix, 'statistics.xlsx'), overwrite = T)

#Counts
write.xlsx(sig_export, file = paste0('Public_data/MCSA_blood/results/', file_prefix, 'DEGene_counts_pval05.xlsx'), overwrite = T)
#write.xlsx(sig_export_padj05, file = paste0('results/', file_prefix, 'DEGene_counts_padj05.xlsx'), overwrite = T)
write.xlsx(check_counts, file = paste0('Public_data/MCSA_blood/results/', file_prefix, 'DS_counts.xlsx'), overwrite = T)
# write.xlsx(test_TPM,file = paste0('Public_data/MCSA_blood/', file_prefix, 'TPM.xlsx'), overwrite = T)


##### BARPLOTS #####
# goi <- sig_res %>% head(20) %>% pull('gene')
goi <- c("IL18R1","IL18RAP","IL17RA","IL17RE","CD248","TGFB1")
goi_data <-  check_counts %>% filter(gene_name %in% goi)
goi_data2 <- rename(goi_data, 
                    # 'APOE33_' = contains('APOE3_3'),
                    # 'APOE34_' = contains('APOE3_4')
                    'Control_' = contains('CN'),
                    'MCI_' = contains('MCI')
                    #'AD_' = contains('dementia')
) 

lineWidth = 1
pointSize = 14
for (i in goi) { 
  goi_data3 <- goi_data2 %>% filter(gene_name == i) %>% pivot_longer(
    cols = -(1:1),
    values_to = c("DS_Counts"),
    names_to = c("condition", "replicate"),
    names_sep  = "_")
  
  goi_data3$condition <- factor(goi_data3$condition, levels=unique(goi_data3$condition))
  
  comparisons <- compare_means(
    data = goi_data3,
    formula = DS_Counts ~ condition,
    method = "t.test",
    p.adjust.method = "BH") %>% filter(p<0.05)
  
  comparison_list_sign <- comparisons %>% mutate(comparison_list = map2(group1, group2,c)) %>% pull(comparison_list)
  
  barplot <- ggplot(goi_data3, aes(x = condition, y = DS_Counts, fill = condition)) + 
    #geom_bar(position = 'dodge', stat = 'summary', fun = mean, width = 0.7, colour="black") +
    #geom_errorbar(stat = 'summary', position = 'dodge', width = 0.4) +
    geom_boxplot() + 
    geom_point(aes(x = condition), position = position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width=0.9)) +
    expand_limits(x = 0, y = 0) +
    theme(panel.spacing = unit(1, "lines")) +
    labs(x = NULL, y = c("Normalized Counts ±SE")) +
    ggtitle(paste0(i))+
    theme(
      text = element_text(size = pointSize, colour = "black"),
      rect = element_blank(),
      line = element_line(size = lineWidth, colour = "black"),
      plot.title  = element_text(color="black", size=20, face="bold.italic"),
      axis.title  = element_text(size = pointSize, colour = "black"),
      axis.text.x  = element_text(size = pointSize , colour = "black", angle = 45, hjust = 1),
      axis.text.y  = element_text(size = pointSize , colour = "black"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = pointSize , colour = "black"),
      # legend.key.height = unit(0.1, "cm"),
      # legend.key.width = unit(0.2, "cm"),
      axis.line = element_line(size = lineWidth, colour = "black"),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
    stat_compare_means(comparisons = comparison_list_sign, method = 't.test', label = "p.signif" ) +
    scale_fill_brewer(palette = "Paired") + scale_y_continuous(expand = c(0, 0, .05, 0))
  
  pdf(file = paste0('Public_data/MCSA_blood/plots/barplots/',file_prefix, i,'.pdf'), pointsize = 10, width = 5, height= 10)
  print(barplot)
  dev.off()
}


#### PATHWAY ANALYSIS #####
file_prefix = 'Female_APOE34_data_MCI_vs_Ctrl_'
pathway_data <- read.xlsx(paste0('Public_data/MCSA_blood/results/',file_prefix, 'DEGene_statistics_pval05.xlsx'))

DEgenes <- pathway_data %>% filter(pvalue < 0.05)  %>% pull('gene_name')
DEgenes <- mapIds(org.Hs.eg.db, DEgenes, 'ENTREZID', 'SYMBOL')

DEG_Log2FC <-  pathway_data %>% filter(pvalue < 0.05) %>%  pull('log2FoldChange')
DEG_Log2FC

names(DEG_Log2FC) = DEgenes
DEG_Log2FC <- sort(DEG_Log2FC, decreasing = TRUE)
DEG_Log2FC

#### CURATED pathways #####
pathways.all.curated <- gmtPathways("GSEA_database/c2.all.v2022.1.Hs.entrez.gmt") 

fgseaRes <- fgsea(pathways=pathways.all.curated, stats=DEG_Log2FC)

fgseaResTidy_curated <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  filter(pval<0.05)

write.xlsx(fgseaResTidy_curated, paste0('Public_data/MCSA_blood/results/Pathway_analysis/FGSEA_enrichment/',file_prefix, 'Curated_pathways.xlsx'))



####BARPLOT FGSEA RESULT ####
file_prefix = 'Female_APOE33_data_MCI_vs_Ctrl_'
pathway_data <- read.xlsx(paste0('Public_data/MCSA_blood/results/Pathway_analysis/FGSEA_enrichment/',file_prefix, 'Curated_pathways.xlsx'), sheet =2)

pathway_data$log_pval <- -log10(pathway_data$pval)
barplot_hallmarks <- ggplot(pathway_data, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES), color = 'black') +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  geom_vline(xintercept = 0, color ='black')+
  theme_gsea+
  scale_fill_gradient2(low = 'blue', mid = "white",
                       high = 'red', midpoint = 0, space = "rgb",
                       na.value = "grey50", guide = "colourbar")
barplot_hallmarks

ggsave(paste0('Public_data/MCSA_blood/plots/pathways/FGSEA_pathways/', file_prefix,'Curated_pathways.pdf'),
       barplot_hallmarks,
       dpi = 300)


###### AGING RELATED SIGNATURE #######
file_prefix = 'Male_APOE33_data_'
pathway_data <- read.xlsx(paste0('Public_data/MCSA_blood/results/',file_prefix, 'DEGene_statistics_pval05.xlsx'))

DEgenes <- pathway_data %>% filter(pvalue < 0.05)  %>% pull('gene_name')
DEG_Log2FC <-  pathway_data %>% filter(pvalue < 0.05) %>%  pull('log2FoldChange')

names(DEG_Log2FC) = DEgenes
DEG_Log2FC <- sort(DEG_Log2FC, decreasing = TRUE)
DEG_Log2FC

#GSEA Aging associated gene pathway analysis
pathways.aging_up <- gmtPathways("GSEA_database/Aging/Aging_UP_genesets.v2023.1.Hs.gmt") 
fgseaRes <- fgsea(pathways=pathways.aging_up, stats=DEG_Log2FC)

fgseaResTidy_aging <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  filter(pval<0.05)

write.xlsx(fgseaResTidy_aging, paste0('results/Pathway_analysis/FGSEA_enrichment/',file_export, 'Hallmark_pathways.xlsx'))


#Public Database Marjolein et al. NATURE COM 2015Aging associated gene pathway analysis
pathways.aging_up <- gmtPathways("GSEA_database/Aging/NAT_COM_2015_Aging_UP_genesets.v2023.2.Hs.gmt") 
fgseaRes <- fgsea(pathways=pathways.aging_up, stats=DEG_Log2FC)

fgseaResTidy_aging <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  filter(pval<0.05)

write.xlsx(fgseaResTidy, paste0('results/Pathway_analysis/FGSEA_enrichment/Aging/',file_export, 'Aging_custom_pathways.xlsx'))





######## SCATTER VOLCANO PLOT ####################
volcano <- read.xlsx(paste0('Public_data/MCSA_blood/results/',file_prefix,'statistics.xlsx'))
volcano$log_pval <- -log10(volcano$pvalue)
volcano$log_baseMean <- log10(volcano$baseMean)


sig_data <- volcano %>% filter(pvalue < 0.05)
nosig_data <- volcano %>% filter(pvalue > 0.05)
UP_data <- sig_data %>% filter(log2FoldChange > 0) 
DOWN_data <- sig_data %>% filter(log2FoldChange < 0)

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

###SELECT DATA
select_UP <- read.xlsx(paste0('Public_data/MCSA_blood/results/',file_prefix,'DEGene_statistics_pval05.xlsx'),sheet = 2) %>% pull('UP')
select_DOWN <- read.xlsx(paste0('Public_data/MCSA_blood/results/',file_prefix,'DEGene_statistics_pval05.xlsx'),sheet = 2) %>% pull('DOWN')

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
  labs(title = str_wrap(paste0(file_prefix),60), x = "Log2FC", y = "-log(pvalue)") +
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


ggsave(paste0('Public_data/MCSA_blood/plots/volcano/',file_prefix, '.pdf'),
       plot = scatter_volcano,
       device = NULL,
       path = NULL,
       #scale = 1,
       width = 20,
       height = 20,
       units = c("cm"),
       dpi = 600)

