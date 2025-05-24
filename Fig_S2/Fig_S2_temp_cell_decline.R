# This script makes a plot that compares how cells react to 10 vs 12 degrees


# Setup -------------------------------------------------------------------


# Load required libraries
library(dplyr)
library(ggplot2)
library(car)
library(ggpubr)
library("RColorBrewer")
library(patchwork)
library(here)


# Set file location
here::i_am("Fig_S2/Fig_S2_temp_cell_decline.R")


# Data Prep ---------------------------------------------------------------

# Load 12 deg sample treatment info
sample_treat_12 <- read.csv(here("Fig_1/sample_treatments.csv")) |> 
  mutate(experiment = "12")

# Load 10 deg sample treatment info
sample_treat_10 <- read.csv(here("Fig_S2/sample_treatments_10.csv")) |> 
  mutate(experiment = "10")


# Load and clean RFU data
RFU_df <- read.csv(here('Fig_S2/RFU10_12.csv')) |> 
  left_join(sample_treat_12, 
            by = c("Sample" = "sample_ID")) |> # add treatments
  left_join(sample_treat_10, by = c("Sample" = "sample_ID")) |> 
  mutate(B12 = factor(coalesce(B12.x, B12.y), # clean artifacts from join
                      levels = c("Y", "N"),
                      labels = c("+B12", "-B12")),
         Temperature = factor(coalesce(temp.x, temp.y),
                              levels = c("4", "10", "12")), 
         Experiment = coalesce(experiment.x, experiment.y),
         Hour = factor(Hour, 
                       levels = c("0",
                                  "24",
                                  "48",
                                  "72")),
         Treatment = paste(Temperature, B12, sep = ", ")) |> 
  select(-ends_with(c(".x", ".y")))
  



# Plotting ----------------------------------------------------------------

# Color palettes for 10 deg
cols_10 <- c("#92C5DE", "#F4A582")
cols_12 <- c("#92C5DE", "#B2182B")
cols_all <-c("#92C5DE", "#F4A582", "#B2182B")


# Boxplot for 12 deg 
bp_12 <- ggplot(
  RFU_df |> filter(Experiment == "12"),
  aes(
    x = Hour,
    y = log(cells_avg),
    grouping = Treatment,
    fill = Temperature,
    alpha = B12
  )
) +
  geom_boxplot() +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_fill_manual(values = cols_12) +
  scale_alpha_manual(values = c(1, .25), labels = c("+", "-")) +
  labs(
    y = "" ,
    x = NULL,
    fill = 'Temperature (°C)',
    alpha = expression("B"[12])
  ) +
  guides(alpha = guide_legend(override.aes = list(
    fill = hcl(c(15, 195), 100, 0, alpha = c(0.25, 1)), colour = NA
  ))) +
  ylim(11, 12.5)


# Boxplot for 10 deg
bp_10 <- ggplot(
  RFU_df |> filter(Experiment == "10"),
  aes(
    x = Hour,
    y = log(cells_avg),
    grouping = Treatment,
    fill = Temperature,
    alpha = B12
  )
) +
  geom_boxplot() +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_fill_manual(values = cols_10) +
  scale_alpha_manual(values = c(1, .25), labels = c("+", "-")) +
  labs(
    y = "" ,
    x = NULL,
    fill = 'Temperature (°C)',
    alpha = expression("B"[12])
  ) +
  guides(alpha = guide_legend(override.aes = list(
    fill = hcl(c(15, 195), 100, 0, alpha = c(0.25, 1)), colour = NA
  ))) +
  ylim(11, 12.5)


# Figs last touches -------------------------------------------------------


# Arrange figs
figure <- ggarrange(bp_12, 
                    bp_10, 
                    ncol = 1, 
                    nrow = 2, 
                    common.legend = TRUE,
                    legend ="bottom", 
                    labels = c("a", "b")) 


(bp_12 + bp_10) + plot_layout(guides = 'collect')

# Annotate with y axis label 
annotated_rfu_boxplot <- annotate_figure(figure,
                left = text_grob(bquote('Estimated Log'[10]~'Cells mL'^-1), 
                                 rot = 90, 
                                 size = 20, 
                                 hjust = .75))

# Save fig to folder
ggsave(annotated_rfu_boxplot, 
       file = here("Fig_S2/Fig_S2_temp_cell_decline.pdf"),
        width = 12, 
       height = 9)





