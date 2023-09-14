# References: https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html

## Gene-level differential expression analysis using DESeq2
## !requireNamespace is used to check whether the package is installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BiocGenerics", force = TRUE)
BiocManager::install("DESeq2", force = TRUE)
BiocManager::install("tidyverse", force = TRUE)
BiocManager::install("tximport", force = TRUE)
BiocManager::install("ggplot2", force = TRUE)
BiocManager::install("biomaRt", force = TRUE)

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
library(umap)

##---------------------------------FILE PREPARATION----------------------------------##

## useMart is a function used to connect to the selected BioMart database and dataset (ListEnsemblArchives() to list host sites)
ensembl_hs_mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
## getBM is a function used to retrieve information from the BioMart database
ensembl_df <- getBM(attributes = c("ensembl_transcript_id_version", "ensembl_gene_id_version", "external_gene_name"), mart = ensembl_hs_mart) %>%
  as_tibble() %>%
  mutate(noVersion = sub('\\..*', '', ensembl_transcript_id_version))%>%
  mutate(noVersion_ENSG = sub('\\..*', '', ensembl_gene_id_version))

# Build tx2gene table from ensembl_df
tx2gene <- ensembl_df %>%
  dplyr::select(noVersion, external_gene_name) %>%
  rename(c('Transcript' = "noVersion" , "Gene" = "external_gene_name" )) %>%
  data.frame()

tx2gene <- read.xlsx('tx2gene.xlsx')

## Create a metadata map using your metadata file.
# metadata <- read.xlsx("metadata.xlsx", sheet = 1) %>% remove_missing()
metadata <- read.xlsx("metadata/metadata_clinical_edited_FINAL.xlsx", sheet = 1)

sample_data <- data.frame(position = metadata$Well, apoe = metadata$Genotyping, mmse = metadata$MMSE, clinical = metadata$Clinical, clinical_2 = metadata$Clinical_2, clinical_3 = metadata$Clinical_3, sex = metadata$SEX, age = metadata$Age, qrds = metadata$QRDS)
sample_data$age_scaled <- scale(sample_data$age)


# Outlier filtering  and subsetting(if necessary)
outliers = c('xxx')
sample_data <- sample_data %>%            # Filter sample_data sheet
  filter(!position %in% outliers)


  
#Analyzis parameters
file_prefix = 'xxx'
experiment = xxx
experiment

## List all directories containing data
#*** All the gene quantification files are in the folder "quant_files", make sure metafile sample order matches order of quant files.
all_files <- list.files(".//transcripts_quant", full.names = T, pattern=NULL, all.files=FALSE)
quant_files <- file.path(all_files, "quant.sf")

position_list <- paste0(experiment$position, collapse = '|')
position_list
sample_files <- grep(position_list, quant_files, value = TRUE)
sample_files

## QC check - Order of metafile samples must match quant file order!  Fix metadata file if not!
all(file.exists(quant_files))
sample_files

# Import Transcript quantification files
txi<-tximport(files=sample_files,type="salmon", tx2gene = tx2gene, ignoreTxVersion = T, countsFromAbundance = "lengthScaledTPM")

## Look at the counts and set the counts as a matrix file.
txicounts<-as.matrix(txi$counts)

# Write the counts to file, round them up, and convert to data.frame and remove txicounts file.
count_data <- txi$counts %>%
  round() %>%
  data.frame()
#rm(txicounts)

# count_batch <- experiment$batch
# count_adjusted <- ComBat_seq(txicounts, batch = count_batch, group = NULL)
# tpm_adjusted <- ComBat_seq(txi$abundance, batch = count_batch, group = NULL)
# mode(count_adjusted) <- 'integer'

#-----------------------------------ANALYSIS BEGINS HERE---------------------------------#
## Create DESeq2Dataset object with variables of interest.
dds <- DESeqDataSetFromTximport(txi, colData = experiment, design = ~ xxx + xxx ) 

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
# dds_result <- results(dds_run, contrast = c('xx','x', 'y'), independentFiltering = T, pAdjustMethod = 'BH')
# dds_result <- lfcShrink(dds_run, contrast = c('xx', 'x', 'y'), res = dds_result, type = 'normal')

#dds_result
summary(dds_result, alpha = 0.05)


## View Total number of normalized counts per sample
dds_run <- estimateSizeFactors(dds_run)
DS_norm_counts <- counts(dds_run, normalized = TRUE)
colSums(DS_norm_counts)

# Transform counts for data visualization
vst <- vst(dds_run, blind=TRUE)

# Plot PCA 
# Add nametags
z <- plotPCA(vst, intgroup=c('clinical'),ntop = 200)
#z + geom_label(aes(label = experiment$position))
theme<-theme(aspect.ratio = 1, panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
pdf(file = paste0('plots/PCA/', file_prefix,'', '.pdf'), pointsize = 10)
z + geom_label(aes(label = experiment$position)) +
  theme + scale_color_brewer(palette = 'Set1') 
dev.off()


# Build significant gene table and extract list of sorted DE genes.
sig_res <- dds_result %>%
  data.frame() %>%
  rownames_to_column(var='gene') %>%    # makes a 'gene' column using column1
  as_tibble %>%
  filter(pvalue < 0.05) %>%
  arrange(pvalue)

sorted_DEGenes <- sig_res$gene

# Export significant gene count data
name_list <- c('gene', (paste0(dds_run$position, "_", dds_run$xxx, "_",dds_run$xxx, "_", dds_run$xxx)))
name_list                        

sig_export <- DS_norm_counts %>%
  data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  as_tibble() %>%
  `colnames<-`(name_list) %>%
  semi_join(sig_res) %>%
  dplyr::slice(match(sorted_DEGenes, gene)) 
# Column Names
check_cols <- c(colnames(sig_export))
check_cols

# #Manual Reordering of Columns (if necessary)
reordered_index <- c(
  grep("xx", names(sig_export), ignore.case = T))

sig_export <- sig_export %>%
  dplyr::select(c(1, reordered_index))

sig_export_2 <- sig_export %>%
  dplyr::select(c(1, reordered_index_2))


check_cols <- c(colnames(sig_export))
check_cols


# Normalized counts
check_counts <- DS_norm_counts %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  dplyr::select(c(1, reordered_index)) %>%
  as_tibble %>% 
  `colnames<-`(check_cols) 

# Results
check_res <- dds_result %>%
  as.data.frame() %>%
  rownames_to_column('gene') 

# Output files
#Statistics
write.xlsx(sig_res, file = paste0('results/', file_prefix, 'DEGene_statistics_pval05.xlsx'), overwrite = T)
write.xlsx(check_res, file = paste0('results/', file_prefix, 'statistics.xlsx'), overwrite = T)

#Counts
write.xlsx(sig_export, file = paste0('results/', file_prefix, 'DEGene_counts_pval05.xlsx'), overwrite = T)
write.xlsx(check_counts, file = paste0('results/', file_prefix, 'DS_counts_check.xlsx'), overwrite = T)


#QC Plots
file_prefix <- 'xxx'
sig_export <- read.xlsx('results/xxx')

goi <- c('xxx')
# goi <- sig_res %>% head(20) %>% pull('gene')
goi_data <-  check_counts %>% filter(gene == goi)
goi_data2 <- rename(goi_data, 
                    'xxx' = contains('xxx')
) 

lineWidth = 1
pointSize = 14
for (i in goi) { 
  goi_data3 <- goi_data2 %>% filter(gene == i) %>% pivot_longer(
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
    geom_bar(position = 'dodge', stat = 'summary', fun = mean, width = 0.7, colour="black") +
    geom_errorbar(stat = 'summary', position = 'dodge', width = 0.4) +
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
      axis.line = element_line(size = lineWidth, colour = "black"),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
    stat_compare_means(comparisons = comparison_list_sign, method = 't.test', label = "p.signif" ) +
    scale_fill_brewer(palette = "Paired") + scale_y_continuous(expand = c(0, 0, .05, 0))
  
  pdf(file = paste0('plots/barplots/',file_prefix, i,'.pdf'), pointsize = 10, width = 5, height= 10)
  print(barplot)
  dev.off()
}



