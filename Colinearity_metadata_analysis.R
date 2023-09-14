###General Statics for metadata
library(caTools)
library(quantmod)
library(MASS)
library(corrplot)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(cowplot)
library(RColorBrewer)
library(sva)
library(mltools)
library(data.table)
library(sjPlot)
library(sjmisc)
library(jtools)
library(sjlabelled)

#Investigation of metadata
data <- read.xlsx(paste0('metadata.xlsx'))
data_select <- data %>% dplyr::select(c(Genotyping,SEX,Age, Clinical, MMSE, QRDS))
model_all <- lm(formula = Clinical ~ Genotyping + SEX + Age + MMSE + QRDS, data = data_select)  
summary(model_all)
colinearity_model <- tab_model(model_all)
colinearity_model
summ(model_all)
ggsave(paste0('plots/statistics_colinearity/','Clinical_metadata_colinearity.pdf'),
       colinearity_model,
       dpi = 600)

#Frequency of age histogram metadata
pointSize= 20
histogram <- ggplot(sample_data, aes(x=age)) + 
  geom_histogram(aes(y=..density..), color="black", fill="lightblue", linetype="dashed", binwidth = 4)+
  geom_density(alpha=.2, fill="cyan",adjust = 2) + 
  theme(text = element_text(size = pointSize, colour = "black"),
        rect = element_blank(),
        line = element_line(size = lineWidth, colour = "black"),
        plot.title  = element_text(color="black", size=20),
        axis.title.x = element_text(size = pointSize , colour = "black"),
        axis.title.y = element_text(size = pointSize , colour = "black"),
        axis.text.x  = element_text(size = pointSize , colour = "black"),
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
  scale_y_continuous(expand = c(0, 0, .05, 0)) +
  scale_x_continuous(expand = c(0, 0, 0, 0), limits=c(0,90)) 
histogram
pdf(file = paste0('plots/histogram/','Histogram_age','.pdf'), pointsize = 10, width = 25, height= 5)
print(histogram)
dev.off()

#Correlation Plot
regression_data <- one_hot(as.data.table(data_select))
M <-cor(regression_data)
head(round(M,2))

dev.off()
corrplot(M, type="lower",order="hclust", tl.col="black", tl.srt=45,sig.level = 0.01, insig = "blank")
ggsave(paste0('plots/correlation_analysis/','Clinical_metadata_analysis.pdf'),
       dpi = 600)


#Barplot for colinearity
lineWidth = 2
pointSize = 30

fit1<- lm(apoe~age)
summary(fit1)

comparisons <- compare_means(
    data = All_data,
    formula = age ~ clinical,
    method = "t.test",
    p.adjust.method = "BH") %>% filter(p<0.05)
  
comparison_list_sign <- comparisons %>% mutate(comparison_list = map2(group1, group2,c)) %>% pull(comparison_list)
comparison_list_sign <- list(c('F','M'))

barplot <- ggplot(xxx, aes(x = clinical_2, y = age, fill = clinical_2)) + 
    geom_boxplot(col ='black') +
    #geom_errorbar(stat = 'summary', position = 'dodge', width = 0.4) +
    geom_point(aes(x = clinical_2), position = position_jitterdodge(jitter.width = 0.3, jitter.height=0.4, dodge.width=0.9)) +
    #expand_limits(x = 0, y = 0) +
    theme(panel.spacing = unit(1, "lines")) +
    labs(x = NULL, y = c("Age")) +
    theme(
      text = element_text(size = pointSize, colour = "black"),
      rect = element_blank(),
      line = element_line(size = lineWidth, colour = "black"),
      plot.title  = element_text(color="black", size=1, face="bold.italic"),
      axis.title  = element_text(size = pointSize, colour = "black"),
      axis.text.x  = element_text(size = pointSize, colour = "black", angle = 45, hjust = 1),
      axis.text.y  = element_text(size = pointSize, colour = "black"),
      axis.ticks = element_line(color = 'black'),
      axis.ticks.length=unit(0.2,"inch"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = pointSize, colour = "black"),
      legend.key.height = unit(0.1, "cm"),
      legend.key.width = unit(0.2, "cm"),
      axis.line = element_line(size = lineWidth, colour = "black"),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
    stat_compare_means(comparisons = comparison_list_sign, method = 't.test', label = "p.signif", size = 10) +
    stat_compare_means(label.y = 100, label.x = 1, size = 10,method = "t.test")+
    scale_fill_brewer(palette = "Paired") + scale_y_continuous(expand = c(0, 1, 0.1, 0))
barplot

dev.off()
pdf(file = paste0('plots/barplots/Colinearity_testing.pdf'), pointsize = 10, height = 10, width = 5)
print(barplot)
dev.off()

Converting age variable into quintiles:
  metadata <- read.xlsx("metadata.xlsx", sheet = 1)
sample_data <- data.frame(position = metadata$Well, apoe = metadata$Genotyping, mmse = metadata$MMSE, clinical = metadata$Clinical, clinical_2 = metadata$Clinical_2, clinical_3 = metadata$Clinical_3, sex = metadata$SEX, age = metadata$Age, qrds = metadata$QRDS)
sample_data$age_scaled <- scale(sample_data$age)

##### Dotplot of Age distribution #######
sample_data 
female_data <- sample_data %>% filter(sex == 'F')
male_date <- sample_data %>% filter(sex == 'M')

dotplot_age <- ggplot(female_data, aes(x = age, y = sex)) + 
  geom_text(data=sample_data[sample_data$sex=='F',], aes(color = clinical_2, size =20), label='\u2640',family = "Arial Unicode MS", position = position_jitterdodge(jitter.width = 0.6, jitter.height=0, dodge.width=0))+
  geom_text(data=sample_data[sample_data$sex=='M',], aes(color = clinical_2, size =20), label='\u2642',family = "Arial Unicode MS", position = position_jitterdodge(jitter.width = 0.6, jitter.height=0, dodge.width=0))+
  theme(panel.spacing = unit(1, "lines")) +
  labs(x = c("Age"), y = NULL) +
  theme(text = element_text(size = pointSize, colour = "black"),
        rect = element_blank(),
        line = element_line(size = lineWidth, colour = "black"),
        plot.title  = element_text(color="black", size=1, face="bold.italic"),
        axis.title  = element_text(size = pointSize, colour = "black"),
        axis.text.x  = element_text(size = pointSize, colour = "black", angle = 45, hjust = 1),
        axis.text.y  = element_text(size = pointSize, colour = "black"),
        axis.ticks = element_line(color = 'black'),
        axis.ticks.length=unit(0.2,"inch"),
        axis.line.y = element_line(color ='white'),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = pointSize, colour = "black"),
        axis.line = element_line(size = lineWidth, colour = "black"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
  scale_color_manual(values = c('grey','red','dodgerblue'))+
  scale_y_discrete(expand = c(0, 0.4, 0.4, 0))
dotplot_age

ggsave(paste0('plots/barplots/xxx','.pdf'),
       dotplot_age,
       dpi = 300,
       width = 25,
       height =10,
       units = 'cm',
       device = cairo_pdf)

