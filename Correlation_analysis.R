#Correlation analysis Protein and RNA
#RNA data
RNA_data <- read.xlsx('results/xxx')
Aurora_genes <- c('xxx')
RNA_data_aurora <- RNA_data %>% 
  #rownames_to_column('gene') %>% 
  filter(gene %in% Aurora_genes) %>% 
  column_to_rownames('gene')

RNA_data_aurora_grouped <- data.frame(
  xxx = rowMeans(select(RNA_data_aurora, contains(("xxx")))))

#Protein data 
all_MFI_files_directory <- list.files(".//AURORA/female_analysis/March_6_analysis/MFI", full.names = T, pattern=NULL, all.files=FALSE)
all_percentage_files_directory <- list.files(".//AURORA/female_analysis/March_6_analysis/percentage", full.names = T, pattern=NULL, all.files=FALSE)

all_MFI_files_names <- gsub('.*/', '', all_MFI_files_directory) 
all_MFI_files_names <- gsub('.xlsx','', all_MFI_files_names)
all_MFI_files_names
all_data <- data.frame()

for (k in all_MFI_files_names) {
  MFI_data <- read.xlsx(paste0('.//AURORA/results/merged_5_clusters/',k,'.xlsx'))
  MFI_data <- MFI_data %>% dplyr::rename('ITGAX'= 'Comp-BUV395-A')
  MFI_data <- MFI_data %>% dplyr::rename('FCGR3A'= 'Comp-BUV496-A')
  MFI_data <- MFI_data %>% dplyr::rename('CD38'= 'Comp-BUV661-A')
  MFI_data <- MFI_data %>% dplyr::rename('PDCD1'= 'Comp-BUV737-A')
  MFI_data <- MFI_data %>% dplyr::rename('IL17A'= 'Comp-BV750-A')
  MFI_data <- MFI_data %>% dplyr::rename('CEACAM8'= 'Comp-FITC-A')
  MFI_data <- MFI_data %>% dplyr::rename('IL10'= 'Comp-BB700-A')
  MFI_data <- MFI_data %>% dplyr::rename('APOE'= 'Comp-PE-A')
  MFI_data <- MFI_data %>% dplyr::rename('TGFB1'= 'Comp-PE-CF594-A')
  MFI_data <- MFI_data %>% dplyr::rename('B2M'= 'Comp-APC-A')
  MFI_data <- MFI_data %>% dplyr::rename('CD14'= 'Comp-Alexa_Fluor_700-A')
  MFI_data <- MFI_data %>% dplyr::rename('HLA-DRA'= 'Comp-APC-H7-A')
  MFI_data <- MFI_data %>% dplyr::rename('ITGAM'= 'Comp-Pacific_Blue-A')
  
  MFI_data <- MFI_data %>% filter(Cluster == 'SUM')
  MFI_data$Group <- k
  MFI_data$group_cluster <- paste0(MFI_data$Group,"-",MFI_data$Cluster)
  MFI_data <- MFI_data[,!(names(MFI_data) %in% c('Group','Cluster'))]
  MFI_data <- MFI_data %>% column_to_rownames('group_cluster')
  
  df <- MFI_data
  all_data <- rbind(all_data, df)
}

Aurora_protein_data <- all_data
Aurora_protein_data <- t(Aurora_protein_data)
Aurora_protein_data <- as.data.frame(Aurora_protein_data)

RNA_data_aurora_grouped
colnames_correlation <- colnames(RNA_data_aurora_grouped)
Aurora_protein_data <- Aurora_protein_data %>%
  `colnames<-`(colnames_correlation)

rownames_correlation <- rownames(RNA_data_aurora_grouped)

Aurora_protein_data <- Aurora_protein_data %>% arrange(rownames(Aurora_protein_data),rownames_correlation)

#Correlation
df1 <- RNA_data_aurora_grouped
df2 <- Aurora_protein_data

sapply(1:nrow(df1), function(i) cor(df1[i,], df2[i,]))

sapply(1:ncol(df1), function(i) cor(df1[,i], df2[,i]))

sapply(seq.int(dim(df1)[1]), function(i) cor(df1[i,], df2[i,]))

correlation <- diag(cor(t(df1), t(df2)))
correlation <- as.data.frame(as.list(correlation))
correlation <- correlation %>% `rownames<-`('R_Correlation')

correlation <- as.data.frame(t(correlation))
correlation <- correlation %>% rownames_to_column('Molecule')

pointSize = 20

dotplot_RNA_Protein_correlation <- ggplot(correlation,aes(x = R_Correlation, y = reorder(Molecule,R_Correlation))) + 
  geom_point(size =10, color = 'black') +
  geom_point(aes(color = R_Correlation), size =9) +
  theme(text = element_text(size = pointSize, colour = "black"),
        rect = element_blank(),
        line = element_line(size = lineWidth, colour = "black"),
        plot.title  = element_text(color="black", size=20),
        axis.title.x = element_text(size = pointSize , colour = "black"),
        axis.title.y = element_blank(),
        axis.text.x  = element_text(size = pointSize , colour = "black", vjust=0.1),
        axis.text.y  = element_text(size = pointSize , colour = "black"),
        axis.ticks = element_line(size = lineWidth, colour = "black"),
        axis.line = element_line(size = lineWidth, colour = "black"),
        legend.position = "right",
        legend.title = element_text(size = pointSize , colour = "black"),
        legend.text = element_text(size = pointSize , colour = "black"),
        legend.key.height = unit(1, "cm"),
        legend.key.width = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.margin=unit(.5, "lines"))+
  scale_color_gradient2(limits=c(-0.5,0.5), low="#112bcc", mid="white", high = "orchid3",na.value = "orchid3") +
  labs(title = str_wrap("Correlation Protein and RNA across conditions", 40), x = 'Correlation Score')


dotplot_RNA_Protein_correlation

dev.off()
pdf(file = paste0('plots/correlation/', 'Protein_and_RNA_correlation_purple_blue_col','.pdf'))
print(dotplot_RNA_Protein_correlation)
dev.off()

#Correlation with MMSE Score
RNA_data <- read.xlsx('xxxx')

Selected_genes <- c('xxx')

RNA_data_aurora <- RNA_data %>% 
  #rownames_to_column('gene') %>% 
  filter(gene %in% Selected_genes) %>% 
  column_to_rownames('gene')

RNA_data_aurora <- as.data.frame(t(RNA_data_aurora))

#SAMPLE_data
metadata <- read.xlsx("metadata.xlsx", sheet = 1)
sample_data <- data.frame(position = metadata$Well, apoe = metadata$Genotyping, mmse = metadata$MMSE, clinical = metadata$Clinical, clinical_2 = metadata$Clinical_2, clinical_3 = metadata$Clinical_3, sex = metadata$SEX, age = metadata$Age,qdrs = metadata$QRDS)
name_list <- c((paste0(sample_data$position, "_", sample_data$clinical_2, "_",sample_data$apoe, "_", sample_data$sex)))
sample_data$name_list <- name_list

sample_data <- sample_data %>% arrange(sample_data$name_list)
sample_data <- sample_data %>% filter(sex =='xxx')
RNA_data_aurora <- RNA_data_aurora %>% arrange(rownames(RNA_data_aurora))

#Sanity check to see if all the files align between metadata and 
all(sample_data$name_list == rownames(RNA_data_aurora))
Correlation_all_data <- cbind(sample_data,RNA_data_aurora)
write.xlsx(Correlation_all_data,paste('results/Correlation.xlsx'))

#Set up model for correlation
model_all <- lm(formula = Correlation_all_data$age ~ Correlation_all_data$xxx)  
summary(model_all)

#Correlation Plot
lineWidth = 1
pointSize = 20
gg_corrlation <- ggplot(cor_data, aes(x= CD274, y =age))+
  geom_point(data = cor_data ,size = 5, color = 'orchid') +
  geom_smooth(data = cor_data, method = "gam", se=TRUE, color="orchid3", formula = y ~ x) +
  stat_cor(data = cor_data, color = 'orchid3', label.x.npc = "left", label.y.npc = "top") +
  theme(text = element_text(size = pointSize, colour = "black"),
        rect = element_blank(),
        line = element_line(size = lineWidth, colour = "black"),
        plot.title  = element_text(color="black", size=20),
        axis.title.x = element_text(size = 15 , colour = "black"),
        axis.title.y = element_text(size = 15 , colour = "black"),
        axis.text.x  = element_text(size = 15 , colour = "black", vjust=0.1),
        axis.text.y  = element_text(size = 15 , colour = "black"),
        axis.ticks = element_line(size = lineWidth, colour = "black"),
        axis.line = element_line(size = lineWidth, colour = "black"),
        legend.position = "right",
        legend.title = element_text(size = pointSize , colour = "black"),
        legend.text = element_text(size = pointSize , colour = "black"),
        legend.key.height = unit(1, "cm"),
        legend.key.width = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.margin=unit(.5, "lines"))

gg_corrlation 

dev.off()
pdf(file = paste0('plots/correlation/xxx/', 'xxx','.pdf'))
print(gg_corrlation)
dev.off()
