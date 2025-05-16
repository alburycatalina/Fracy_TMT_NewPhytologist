# This script creates a boxplot with protein quotas for figure S3 

# Load required pacakges
library(tidyverse)
library(here)

# Load in data 
prot_quota_df <- read.csv(here("Fig_2/Fig_2_protein_quota_data/prot_quota_df.csv")) |> 
  
  # set treatments as levels
  mutate(Temperature = factor(Temperature, levels = c("4", "12")),
         B12 = factor(B12, 
                      levels = c("+B12", "-B12")))

# Protein Boxplot (Fig S3) ---------------------------------------------------------

# Color palette 
col_12 <- c("#92C5DE", "#B2182B" )

# Boxplot
protquota_bp1 <-
  ggplot(
    prot_quota_df,
    aes(
      x = Temperature,
      y = pg_prot_cell,
      grouping = Treatment,
      fill = Temperature,
      alpha = B12 # Make +B12 more saturated than -B12 in color
    )
  ) +
  geom_boxplot() +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_fill_manual(values = col_12, labels = c("4", "12")) +
  scale_alpha_manual(values = c(1, .25), labels = c("+", "-")) +
  labs(
    y = bquote('Total Protein (Picograms cell' ^ -1 * ')') ,
    x = "Treatment",
    fill = 'Temperature (°C)',
    alpha = expression("B"[12])
  ) +
  guides(alpha = guide_legend(override.aes = list(
    fill = hcl(c(15, 195), 100, 0, alpha = c(0.25, 1)), colour = NA
  ))) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())


# Save plot as pdf 
ggsave(here("Fig_S3/Fig_S3_total_prot_boxplot.pdf"), 
       width = 7, 
       height = 6, 
       units = "in")
