
# This script produces a plot for cell counts and fvfm over time after exposure to temperature treatment. 

library(tidyverse)
library(here)
library(broom)
library(tibble)
library(ggpubr)


# Cell Count Data Cleaning ------------------------------------------------

here::i_am("Fig_01/Fracy_TMT_Fig_01_cellcnts_fvfm.R")


# Clean cell cnt data
cell_cnts <- read.csv(here("Fig_01/Fig_01_RAW/cell_cnt_raw.csv")) |> 
  
  # exclude 96 hour data 
  # lack of new media makes 4deg start to decline
  filter(time_since_temp_hours != "96") |> 
  
  # convert some data to factors
  mutate(time_since_temp_hours = factor(time_since_temp_hours),
        Temperature = factor(Temperature), 
        B12 = factor(B12))



# Cell Count Barplot --------------------------------------------------------------------
# Set colors
temp_palette <- c("#92C5DE", "#B2182B")

# Barplot
cell_cnt_barplot <- ggplot(
  cell_cnts,
  aes(
    x = time_since_temp_hours,
    y = cells_mL,
    grouping = Treatment,
    fill = Temperature,
    alpha = B12
  )
) +
  geom_boxplot() +
  theme_classic() +
  theme(text = element_text(size = 22)) +
  scale_fill_manual(values = temp_palette) +
  scale_alpha_manual(values = c(1, .25), labels = c("+", "-")) +
  labs(
    y = bquote('Cells mL' ^ -1),
    x = 'Time (Hours)',
    fill = 'Temperature (°C)',
    alpha = expression("B"[12])
  ) +
  guides(alpha = guide_legend(override.aes = list(
    fill = hcl(c(15, 195), 100, 0, alpha = c(0.25, 1)), colour = NA
  ))) +
  theme(
    # panel.background = element_rect(fill = "transparent"),
    # # bg of the panel
    # plot.background = element_rect(fill = "transparent", color = NA),
    # # bg of the plot
    # panel.grid.major = element_blank(),
    # # get rid of major grid
    # panel.grid.minor = element_blank(),
    # # get rid of minor grid
    # legend.background = element_rect(fill = "transparent"),
    # # get rid of legend bg
    # legend.box.background = element_rect(fill = "transparent"),
    # legend.key = element_blank(),
    # # get rid of legend panel bg,
    text = element_text(size = 25)
  ) 


# Student's T-test Cell Counts -------------------------------------------------------------------

# Function for getting t-test against different time points 

cell_cnt_ttest <- function(df, hour){
  
  # x
  x <- df |> 
    filter(time_since_temp_hours == paste(hour) & 
             Temperature == "4") |> 
    pull(cells_mL)
  
  
  y <- df |>  
    filter(time_since_temp_hours == paste(hour) & 
             Temperature == "12") |> 
    pull(cells_mL)
  
  assign(paste("ttest", hour, sep = "_"),
         
         # Complete t-test  
         tidy(t.test(x, y)), envir = globalenv())
          
  
}

# Get list of time points in cell count data 
hours <- levels(cell_cnts$time_since_temp_hours)

# Apply function to get t-test stats for all hours 
cellcnts_stats_output <- sapply(hours, cell_cnt_ttest, df = cell_cnts)  |>
    t() |> 
  as.data.frame() |> 
  rownames_to_column(var = "time_since_temp_hours")
  


# Clean up fv/fm data -----------------------------------------------------

# Read in table with sample treatment info 
sample_treat <- read.csv(here("Fig_01/Fig_01_RAW/sample_treatments.csv"))

# Read in and clean fvfm data 
fvfm <- read.csv(here("Fig_01/Fig_01_RAW/fvfm_raw.csv")) |>
  rename(c( # fix column headers
    "sample_ID" = "Sample",
    "Hour" = "h",
    "FvFm_1" = "FvFm",
    "FvFm_2" = "X",
    "fvfm_mean" = "AVG"
  )) |> 
  left_join(sample_treat, by = "sample_ID") |> 
  mutate(temp = factor(temp),
         Hour = factor(Hour))


# Fv/Fm Barplot -----------------------------------------------------------


fvfm_barplot <- ggplot(fvfm,
                       aes(
                         x = Hour,
                         y = fvfm_mean,
                         grouping = interaction(temp, B12),
                         fill = temp,
                         alpha = B12
                       )) +
  geom_boxplot() +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_fill_manual(values = temp_palette) +
  scale_alpha_manual(values = c(1, .25), labels = c("+", "-")) +
  labs(
    y = bquote('F'[v] ~ '/F'[m]),
    x = 'Time (Hours)',
    fill = 'Temperature (°C)',
    alpha = expression("B"[12])
  ) +
  guides(alpha = guide_legend(override.aes = list(
    fill = hcl(c(15, 195), 100, 0, alpha = c(0.25, 1)), colour = NA
  ))) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
  theme(
    # panel.background = element_rect(fill = "transparent"),
    # # bg of the panel
    # plot.background = element_rect(fill = "transparent", color = NA),
    # # bg of the plot
    # panel.grid.major = element_blank(),
    # # get rid of major grid
    # panel.grid.minor = element_blank(),
    # # get rid of minor grid
    # legend.background = element_rect(fill = "transparent"),
    # # get rid of legend bg
    # legend.box.background = element_rect(fill = "transparent"),
    # get rid of legend panel bg,
    text = element_text(size = 25)
  )


# Arrange and Export Plots -----------------------------------------------------------
# Arrange plots
barplots <- ggarrange(
  cell_cnt_barplot,
  fvfm_barplot,
  ncol = 2,
  nrow = 1,
  common.legend = TRUE, 
  # individual plot labels
  labels = c("a", "b"),  
  font.label = list(size = 25, 
                    color = "black"),
  legend = "bottom"
)

# Export plot
ggsave(barplots, 
       filename = here("Fig_01/Fig_01_cellcnts_fvfm_barplot.pdf"),  
       bg = "transparent", 
       width = 15,
       height = 9)
