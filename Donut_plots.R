#FACS DONUT PLOTS
library(ggplot2)
library(dplyr)

#Download data 
FACS_data <- read.xlsx('results/xxx', sheet =2)

FACS_data_xxx <- FACS_data %>% filter(Group == 'xxx')
FACS_data_xx <- FACS_data %>% filter(Group == 'xx')

label <- 'xxx'
experiment_FACS <- 'xxx'

# Hole size
hsize <- 3

color_clusters <- c('mediumpurple','cyan3','orangered','plum1')

donut <- ggplot(xxx, aes(ymax=10, ymin=8, xmax=4, xmin=3, x = hsize, y = Percentage, fill = Label)) +
  geom_col(color = "black") +
  geom_text(aes(label = Percentage),
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = color_clusters) +
  xlim(c(0.2, hsize + 0.5)) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  ggtitle(paste0(label," ",experiment_FACS))
donut

color_clusters <- c('mediumpurple','cyan3','orangered','plum1')

hsize_2 <- 2
donut_2 <- ggplot(xxx, aes(x = hsize_2, y = Percentage, fill = Label)) +
  geom_col(color = "black") +
  geom_text(aes(label = Percentage),
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = color_clusters) +
  xlim(c(0.2, hsize + 0.5)) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  ggtitle(paste0(label," ",experiment_FACS))
donut_2

ggsave(paste0('plots/FACS_donutplots/',label,'_', experiment_FACS,'.pdf'),
         plot = donut,
         device = NULL,
         path = NULL,
         #scale = 1,
         width = 20,
         height = 20,
         units = c("cm"),
         dpi = 300)
  
