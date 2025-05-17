# This script produces figure 2, which contains information on the amount of protein per cell 



# Setup -------------------------------------------------------------------

# Load required packages
library(ggplot2)
library(agricolae)
library(dplyr)
library(ggpubr)
library(here)
library(forcats)

# Declare file locations with `here`
here::i_am("Fig_2/Fracy_TMT_Fig_2_methionine_synthase_SRM.R")



# Protein Quota Calculations ----------------------------------------------------------

# Load in cell counts
cell_cnt_data <- read.csv(here("Fig_1/cell_cnt_raw.csv")) |> 
  filter(time_since_temp_hours == 24) |> # hour 24 is time of harvest 
  select(Sample, B12, Temperature, Treatment, cells_mL) |>
  mutate(Sample = substring(Sample, 5)) |>   # remove cell numbers 
  rename(sample_id = Sample)

# Get protein quotas 
# Use BCA run to find the concentration of the protein before digestion for MS 
prot_quota_df <- 
  read.csv(here(
  "Fig_2/Fig_2_Raw_data/BCA1.1_05062019_a.txt_dilution25.csv")) |> # load in protein BCA data 
  select(sample_id, ug_extracted_total) |>
  mutate(
    sample_id = c(
      'CMA-1-119-1',
      'CMA-1-119-2',
      'CMA-1-119-3',
      'CMA-1-137-7',
      'CMA-1-127-8',
      'CMA-1-119-9',
      'CMA-1-146-1',
      'CMA-1-146-2',
      'CMA-1-146-3',
      'CMA-1-146-7',
      'CMA-1-146-8',
      'CMA-1-146-9'
    )
  ) |> 
  full_join(cell_cnt_data, by = "sample_id") |> 
  mutate(cells_filter = cells_mL* 20, # 20 mL of culture filtered (2 x 10 mL filters in this extraction)
         ugTotalProt_cell = ug_extracted_total / cells_filter, # calculate protein per cell
         pgTotalProt_cell = ugTotalProt_cell * 10^6
  ) |> 
  
  # set treatments as levels
  mutate(Temperature = factor(Temperature, levels = c("4", "12")),
         B12 = factor(B12, 
                      levels = c("Y", "N"),
                      labels = c("+B12", "-B12")))

# Write protein quota data 
write.csv(prot_quota_df, file = here("Fig_2/Fig_2_protein_quota_data/prot_quota_df.csv"), 
          row.names = FALSE)


# Protein SRM Calculations ------------------------------------------------

# Calculate pmol/ug total protein and molecules per cell of analytes from SRM
targeted_data_raw <- read.csv(here('Fig_2/Fig_2.csv')) |> # Load in protein data
  
  rename(sample_id = Harvest_ID) |> 
  
  # Calculate ratio of light to heavy peptide
  mutate(light_heavy_ratio = peak_area_light/peak_area_heavy) |> 
  
  # Calculate fmol of each protein per total ug of protein from digest (20 ug of protein was digested)
  mutate(fmolAnalyte_ugProtein = light_heavy_ratio * 20) |> 
  
  # Calculate fmol of protein per ug protein on column
  mutate(fmolAnalyte_ugProtein = fmolAnalyte_ugProtein/ug_on_col,
         pmolAnalyte_ugProtein = fmolAnalyte_ugProtein / 10^3) |> 
  
  # Bring in protein quota data from 
  left_join({prot_quota_df |> 
              select(sample_id, ugTotalProt_cell)}, 
            by = "sample_id") |> 
  
  # Calculate molecules per cell 
  mutate(pmolAnalyte_cell = pmolAnalyte_ugProtein * ugTotalProt_cell, 
         molAnalyte_cell = pmolAnalyte_cell * 10^15, 
         
         # Multiply by avagadro's number for molecules
         moleculesAnalyte_cell = molAnalyte_cell * 6.022E23)
  
  

# Summary Stats -----------------------------------------------------------


# Roll up to mean pmol per ug protein for replicates
targeted_data_raw_repsum <- targeted_data_raw |> 
  group_by(sample_id, 
           B12, 
           Temperature, 
           Protein) |> 
  dplyr::summarise(repmean_pmolAnalyte_ugProtein = mean(pmolAnalyte_ugProtein),
                   repsd_pmolAnalyte_ugProtein = sd(pmolAnalyte_ugProtein))

# Roll up to mean pmol per ug protein for treatments
targeted_dataraw_treatsum <- targeted_data_raw_repsum |> 
  group_by(B12, 
           Temperature, 
           Protein) |> 
  dplyr::summarise(treatmean_pmolAnalyte_ugProtein = mean(repmean_pmolAnalyte_ugProtein),
                   treatsd_pmolAnalyte_ugProtein = sd(repmean_pmolAnalyte_ugProtein))

write.csv(targeted_dataraw_treatsum, 
          file = here("Fig_2/Fig_2_protein_quota_data/targeted_dataraw_treatsum.csv"))




# Data Prep for plots ---------------------------------------------------------------

# Create new df for plots
targeted_data <- targeted_data_raw |>
  
  # Begin to convert to long format
  dplyr::select(sample_id,
                Protein,
                B12,
                Temperature,
                pmolAnalyte_ugProtein,
                moleculesAnalyte_cell) |>
  
  # Convert to a wide format for plotting
  pivot_wider(
    values_from = c(pmolAnalyte_ugProtein, moleculesAnalyte_cell),
    names_from = Protein,
    names_glue = "{Protein}_{.value}"
  ) |>
  
  # Fold change calculations
  mutate(
    fc = (MetE_pmolAnalyte_ugProtein - MetH_pmolAnalyte_ugProtein) / MetH_pmolAnalyte_ugProtein
  ) |>
  
  # Set up levels for plotting
  mutate(
    Temperature = factor(Temperature, 
                         levels = c("4", "12")),
    B12 = factor(paste(B12, "B12", sep = ""), 
                 levels = c("+B12", "-B12")),
    Treatment = factor(paste0(Temperature, ",", B12),
                       levels = c("4,+B12", 
                                  "4,-B12", 
                                  "12,+B12", 
                                  "12,-B12")
  ))


# MetH/MetE Boxplots -------------------------------------------------------------------

# Color palette 
col_12 <- c("#92C5DE", "#B2182B" )

# MetH Plot
metH_quota_picomol <- ggplot(
  targeted_data,
  aes(
    x = Temperature,
    y = MetH_pmolAnalyte_ugProtein,
    grouping = Treatment,
    fill = Temperature,
    alpha = B12 # Make B12 dictate the transparency 
  )
) +
  geom_boxplot() + # Initiate a boxplot
  theme_classic() + # Classic theme
  
  # Make text larger
  theme(text = element_text(size = 20)) +
  
  # Set colors from palette
  scale_fill_manual(values = col_12, 
                    labels = c("4",
                               "12")) +
  
  # Make B12 the transparency 
  scale_alpha_manual(
    values = c(1, .25),
    labels = c("+", "-"),
    breaks = c("+B12", "-B12")
  ) +
  
  # Labels
  labs(
    y = bquote('Picomoles MetH \u03bcg protein' ^ -1) , # ilu unicode characters <3 
    x = "Treatment",
    fill = 'Temperature (°C)',
    alpha = expression("B"[12])
  ) +
  
  # Fix legend
  guides(alpha = guide_legend(override.aes = list(
    fill = hcl(c(15, 195), 
               100,
               0, 
               alpha = c(.25, 1)), 
    colour = NA
  ))) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  ylim(0, .025)

# MetE Plot
metE_quota_picomol <- ggplot(
  targeted_data,
  aes(
    x = Temperature,
    y = MetE_pmolAnalyte_ugProtein,
    grouping = Treatment,
    fill = Temperature,
    alpha = B12
  )
) +
  geom_boxplot() +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_fill_manual(values = col_12, labels = c("4", "12")) +
  scale_alpha_manual(values = c(1, .25), labels = c("+", "-")) +
  labs(
    y = bquote('Picomoles MetE \u03bcg protein' ^ -1) ,
    x = "Treatment",
    fill = 'Temperature (°C)',
    alpha = expression("B"[12])
  ) +
  guides(alpha = guide_legend(override.aes = list(
    fill = hcl(c(15, 195), 100, 0, alpha = c(0.25, 1)), colour = NA
  ))) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  ylim(0, .025)


# MetH Plot
metH_quota_molecules <- ggplot(
  targeted_data,
  aes(
    x = Temperature,
    y = MetH_moleculesAnalyte_cell,
    grouping = Treatment,
    fill = Temperature,
    alpha = B12 # Make B12 dictate the transparency 
  )
) +
  geom_boxplot() + # Initiate a boxplot
  theme_classic() + # Classic theme
  
  # Make text larger
  theme(text = element_text(size = 20)) +
  
  # Set colors from palette
  scale_fill_manual(values = col_12, 
                    labels = c("4",
                               "12")) +
  
  # Make B12 the transparency 
  scale_alpha_manual(
    values = c(1, .25),
    labels = c("+", "-"),
    breaks = c("+B12", "-B12")
  ) +
  
  # Labels
  labs(
    y = bquote('Molecules MetH cell' ^ -1) , # ilu unicode characters <3 
    x = "Treatment",
    fill = 'Temperature (°C)',
    alpha = expression("B"[12])
  ) +
  
  # Fix legend
  guides(alpha = guide_legend(override.aes = list(
    fill = hcl(c(15, 195), 
               100,
               0, 
               alpha = c(.25, 1)), 
    colour = NA
  ))) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  ylim(0, 3e32)

# MetE Plot
metE_quota_molecules <- ggplot(
  targeted_data,
  aes(
    x = Temperature,
    y = MetE_moleculesAnalyte_cell,
    grouping = Treatment,
    fill = Temperature,
    alpha = B12
  )
) +
  geom_boxplot() +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_fill_manual(values = col_12, labels = c("4", "12")) +
  scale_alpha_manual(values = c(1, .25), labels = c("+", "-")) +
  labs(
    y = bquote('Molecules MetE cell' ^ -1) ,
    x = "Treatment",
    fill = 'Temperature (°C)',
    alpha = expression("B"[12])
  ) +
  guides(alpha = guide_legend(override.aes = list(
    fill = hcl(c(15, 195), 100, 0, alpha = c(0.25, 1)), colour = NA
  ))) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  ylim(0, 3e32)


# Arrange MetH/MetE Boxplots
ggarrange(metH_quota_picomol, 
          metE_quota_picomol,
          metH_quota_molecules,
          metE_quota_molecules,
          ncol = 2,
          nrow = 2, 
          common.legend = TRUE, 
          legend = "bottom", 
          labels = c("a", "b", "c", "d"),  
          font.label = list(size = 25, 
                            color = "black"))

# No 
ggsave(here("Fig_2/Fig_2_meth_mete_barplots.png"), 
       width = 12,
       height = 12, 
       units = "in")

# t-test on MetH ------------------------------------------------------

# Temp MetH - not sig (p = 0.08816)
t.test({targeted_data |> 
         filter(Temperature == "4") |> 
         pull(MetH_pmolAnalyte_ugProtein)}, 
       {targeted_data |> 
         filter(Temperature == "12") |> 
         pull(MetH_pmolAnalyte_ugProtein)})


# B12 MetH - not sig (p-value = 0.615)
t.test({targeted_data |> 
    filter(B12 == "+B12") |> 
    pull(MetH_pmolAnalyte_ugProtein)}, 
    {targeted_data |> 
        filter(B12 == "-B12") |> 
        pull(MetH_pmolAnalyte_ugProtein)})


# t-test on MetE ------------------------------------------------------

# Temp MetE - not sig (p-value = 0.5993)
t.test({targeted_data |> 
    filter(Temperature == "4") |> 
    pull(MetE_pmolAnalyte_ugProtein)}, 
    {targeted_data |> 
        filter(Temperature == "12") |> 
        pull(MetE_pmolAnalyte_ugProtein)})



# B12 MetE - significantly less in noB12 samples (p-value = 0.0003055)
t.test({targeted_data |> 
    filter(B12 == "+B12") |> 
    pull(MetE_pmolAnalyte_ugProtein)}, 
    {targeted_data |> 
        filter(B12 == "-B12") |> 
        pull(MetE_pmolAnalyte_ugProtein)})



