# This script plots a curve culture and RFU x cell count relationships

# Load required packages
library(here)
library(tidyverse)
library(ggpubr)


# Load in data fror curve culture
curve_data <- read.csv('curve_culture.csv')

curve_1 <- ggplot(curve_data, aes(x= Day, y=RFU)) +
  geom_path() +
  xlab('Time (Days)') +
  ylab('Fluorescence (RFU)') +
  theme_classic() +
  theme(axis.text=element_text(size=15, color = 'black'), 
        axis.title = element_text(size= 17, face='bold')) +
  geom_errorbar(aes(ymin= RFU-stdev, ymax= RFU+stdev), width= .4)

curve_2 <- ggplot(curve_data, aes(x=Day, y=log(RFU))) +
  geom_path() +
  xlab('Time (Days)') +
  ylab('ln(RFU)') +
  theme_classic() +
  geom_vline(xintercept = 5, linetype="dashed", color = "red") +
  geom_vline(xintercept = 10, linetype="dashed", color = "red") +
  annotate('text', x= 7.5, y=5.5, label= 'Logistical Growth') +
  annotate('text', x= 7.5, y=5.35, label= 'R² = 0.9998') +
  theme(axis.text=element_text(size=15, color = 'black'), 
        axis.title = element_text(size= 17, face='bold'))

log_phase <- curve_data[1:3,]

summary(lm(Day~log(RFU), data= log_phase))



# Load in data from curve culture
curve <- read.csv('RFU_cellcnt_curve.csv') |> 
  mutate(Mean_RFU = as.numeric(Mean_RFU))

curve_line <- lm(curve$Mean_RFU ~ curve$Cell_Count) 
summary(curve_line)

cells_rfu_lm <- ggplot(curve, 
       aes(x = Mean_RFU, 
           y = Cell_Count)) +
  geom_point(size = 5) +
  geom_smooth(color = 'black', 
              method = 'lm', 
              se = FALSE) +
  theme_classic() +
  ylab(bquote('Cells mL'^-1)) +
  xlab('RFU')+
  theme(axis.text=element_text(size=15, color = 'black'), 
        axis.title = element_text(size= 20))
  
# Facet plots together 
curve_plots <- ggarrange(curve_1, 
          curve_2, 
          cells_rfu_lm, 
          nrow = 1, 
          ncol = 3, 
          labels = c("a", "b", "c"),
          font.label = list(size = 25, 
                            color = "black"))

ggsave(curve_plots, 
       file = here("Fig_S1/Fig_S1_curve_cultures.pdf"), 
       width = 20, 
       height = 8)

