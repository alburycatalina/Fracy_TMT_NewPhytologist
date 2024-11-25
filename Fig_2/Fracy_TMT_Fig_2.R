# This script produces figure 2, which contains information on the amount of protein per cell 



# Setup -------------------------------------------------------------------
setwd("~/Bertrand Lab Dropbox/Bertrand Lab shared workspace/Catalina/Summer_2022/0_Fracy_TMT_Manuscript/Fracy_TMT_FigScripts/Fig_2")

# Load required packages
library(ggplot2)
library(agricolae)
library(dplyr)
library(ggpubr)


# Protein SRM Calculations ------------------------------------------------


# Load in protein data
targeted_data_raw <- read.csv('Fig_2.csv') |> 
  # Calculate ratio of light to heavy peptide
  mutate(light_heavy_ratio = peak_area_light/peak_area_heavy) |> 
  
  # Calculate fmol of each protein per total ug of protein from digest (20 ug of protein was digested)
  mutate(fmol_ug_prot = light_heavy_ratio*20) |> 
  
  # Calculate fmol of protein per ug protein on column
  mutate(fmolProtein_ugProtein = fmol_ug_prot/ug_on_col) |> 
  
  # Calculate fmol protein per cell 
  mutate(fmolProtein_Cell = fmolProtein_ugProtein * ug_total_prot_cell) 




# Data Prep for plots ---------------------------------------------------------------

# Create new df for plots
targeted_data <- targeted_data_raw |> 
  
  # Begin to convert to long format 
  dplyr::select(Harvest_ID, Protein, B12, Temperature, ug_total_prot_cell) |> 
  filter(Protein == "MetH") |> 
  select(-c(Protein)) |> 
  mutate(MetH_fmol_cell = targeted_data_raw$fmolProtein_Cell[targeted_data_raw$Protein == "MetH"], 
         MetE_fmol_cell = targeted_data_raw$fmolProtein_Cell[targeted_data_raw$Protein == "MetE"],
         ACT1_fmol_cell = targeted_data_raw$fmolProtein_Cell[targeted_data_raw$Protein == "ACT1"],
         RBCL_fmol_cell = targeted_data_raw$fmolProtein_Cell[targeted_data_raw$Protein == "RBCL"]) |> 
  
  # Convert to yoctomoles
  mutate(MetH_ymol_cell = MetH_fmol_cell * 1e9,
         MetE_ymol_cell = MetE_fmol_cell * 1e9,
         ACT1_ymol_cell = ACT1_fmol_cell * 1e9,
         RBCL_ymol_cell = RBCL_fmol_cell * 1e9,
         ACT1_ymol_cell = ACT1_fmol_cell * 1e9,
         RBCL_ymol_cell = RBCL_fmol_cell * 1e9) |>  
  
  
  # Fold change calculations
  mutate(fc = (MetE_fmol_cell - MetH_fmol_cell)/MetH_fmol_cell) |> 

  # Set up levels for plotting
  mutate(Temperature = factor(Temperature,
                              levels = c("4", 
                                         "12")),
         B12 = factor(paste(B12, 
                            "B12", 
                            sep= ""),
                      levels = c("+B12", 
                                 "-B12")),
         Treatment = factor(paste(Temperature,
                                  ",",
                                  B12)))


# Summary Stats -----------------------------------------------------------


# Roll up to mean fmol per ug protein
targeted_data_raw_repsum <- targeted_data_raw |> 
  group_by(Harvest_ID, 
           B12, 
           Temperature, 
           Protein) |> 
  dplyr::summarise(mean_fmol_ug = mean(fmolProtein_ugProtein))

targeted_dataraw_treatsum <- targeted_data_raw_repsum |> 
  group_by(B12, 
           Temperature, 
           Protein) |> 
  dplyr::summarise(mol_ug = mean(mean_fmol_ug),
                   sd_fmol_ug = sd(mean_fmol_ug))
  
  
# MetH/MetE Boxplots -------------------------------------------------------------------

# Color palette 
col_12 <- c("#92C5DE", "#B2182B" )

# MetH Plot
metH_quota <- ggplot(
  targeted_data,
  aes(
    x = Temperature,
    y = MetH_ymol_cell,
    grouping = Treatment,
    fill = Temperature,
    alpha = B12 # Make B12 the transparency 
  )
) +
  geom_boxplot() + # Initate a boxplot
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
    y = bquote('Yoctomoles MetH cell' ^ -1) ,
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
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())

# MetE Plot
metE_quota_bg <- ggplot(
  targeted_data,
  aes(
    x = Temperature,
    y = MetE_ymol_cell,
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
    y = bquote('Yoctomoles MetE cell' ^ -1) ,
    x = "Treatment",
    fill = 'Temperature (°C)',
    alpha = expression("B"[12])
  ) +
  guides(alpha = guide_legend(override.aes = list(
    fill = hcl(c(15, 195), 100, 0, alpha = c(0.25, 1)), colour = NA
  ))) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) 

# Arrange MetH/MetE Boxplots
ggarrange(metH_quota, 
          metE_quota_bg, 
          ncol = 3,
          nrow = 1, 
          common.legend = TRUE, 
          legend = "bottom", 
          labels = c("a", "b"),  
          font.label = list(size = 25, 
                            color = "black"))

# Transparent Background Figs ---------------------------------------------
# MetH Plot
metH_quota_bg <- ggplot(targeted_data, aes(x=Temperature, 
                                           y= MetH_fmol_cell, 
                                           grouping = Treatment, 
                                           fill = Temperature, 
                                           alpha = B12)) + 
  geom_boxplot() +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_fill_manual(values = col_12, labels= c("4", "12")) +
  scale_alpha_manual(values = c(1, .25), 
                     labels = c( "+", "-")) +
  labs(y= bquote('Femtomoles MetH cell'^-1) , 
       x= "Treatment", 
       fill='Temperature (°C)', 
       alpha = expression("B"[12])) +
  guides(alpha=guide_legend(override.aes=list(fill=hcl(c(15,195),100,0, alpha=c(0.25,1)), colour=NA))) +
  theme(axis.text.x = element_blank(), 
        axis.ticks = element_blank()) +
  ylim(0, 8e-9) +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg, 
    text = element_text(size=25)
  )


ggsave(metH_quota_bg, filename = "metH_transbg.png",  bg = "transparent")

metE_quota_bg <- ggplot(
  targeted_data,
  aes(
    x = Temperature,
    y = MetE_fmol_cell,
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
    y = bquote('Femtomoles MetE cell' ^ -1) ,
    x = "Treatment",
    fill = 'Temperature (°C)',
    alpha = expression("B"[12])
  ) +
  guides(alpha = guide_legend(override.aes = list(
    fill = hcl(c(15, 195), 100, 0, alpha = c(0.25, 1)), colour = NA
  ))) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  theme(
    panel.background = element_rect(fill = "transparent"),
    # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA),
    # bg of the plot
    panel.grid.major = element_blank(),
    # get rid of major grid
    panel.grid.minor = element_blank(),
    # get rid of minor grid
    legend.background = element_rect(fill = "transparent"),
    # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"),
    # get rid of legend panel bg,
    text = element_text(size = 25)
  )

ggsave(metE_quota_bg, filename = "metE_transbg.png",  bg = "transparent")


# t-test on MetH ------------------------------------------------------

# Temp MetH - not sig (p = 0.6726)
x_methTemp <- targeted_data |> 
  filter(Temperature == "4") |> 
  pull(MetH_fmol_cell)

y_methTemp <- targeted_data |> 
  filter(Temperature == "12") |> 
  pull(MetH_fmol_cell)

t.test(x_methTemp, y_methTemp)


# B12 MetH - not sig (p-value = 0.5063)
x_methB12 <- targeted_data |> 
  filter(B12 == "+B12") |> 
  pull(MetH_fmol_cell)

y_methB12 <- targeted_data |> 
  filter(B12 == "-B12") |> 
  pull(MetH_fmol_cell)

t.test(x_methB12, y_methB12)


# t-test on MetE ------------------------------------------------------

# Temp MetE - not sig (p-value = 0.6615)
x_meteTemp <- targeted_data |> 
  filter(Temperature == "4" & 
           B12 == "-B12") |> 
  pull(MetE_fmol_cell)

y_meteTemp <- targeted_data |> 
  filter(Temperature == "12" & 
           B12 == "-B12") |> 
  pull(MetE_fmol_cell)

t.test(x_meteTemp, y_meteTemp)


# B12 MetE - significantly less in noB12 samples (p-value = 0.001651)
x_meteB12 <- targeted_data |> 
  filter(B12 == "+B12") |> 
  pull(MetE_fmol_cell)

y_meteB12 <- targeted_data |>  
  filter(B12 == "-B12") |>
  pull(MetE_fmol_cell)

t.test(x_meteB12, y_meteB12)

# Break down by temps
x_meteB12_4 <- targeted_data |> 
  dplyr::filter(B12 == "+B12" &
                  Temperature == "4") |> 
  pull(MetE_fmol_cell)

y_meteB12_4 <- targeted_data |>  
  filter(B12 == "-B12" & 
           Temperature == "4") |> 
  pull(MetE_fmol_cell)

t.test(x_meteB12_4, y_meteB12_4)

# 12
x_meteB12_12 <- targeted_data |> 
  dplyr::filter(B12 == "+B12" &
                  Temperature == "12") |> 
  pull(MetE_fmol_cell)

y_meteB12_12 <- targeted_data |>
  filter(B12 == "-B12" & 
           Temperature == "12") |> 
  pull(MetE_fmol_cell)


t.test(x_meteB12_12, y_meteB12_12)

# Protein Quotas ----------------------------------------------------------

# load in protein quota data
protquota_data <- read.csv("Fig_2_Raw_data/BCA1.1_05062019_a.txt_dilution25.csv")

#enter 12deg samples treatment data
id_12 <- c('CMA-1-119-1','CMA-1-119-2','CMA-1-119-3', 'CMA-1-137-7', 'CMA-1-129-8', 'CMA-1-119-9', 'CMA-1-146-1', 'CMA-1-146-2', 'CMA-1-146-3', 'CMA-1-146-7', 'CMA-1-146-8', 'CMA-1-146-9')

B12_12 <- c('+B12','+B12','+B12','-B12','-B12','-B12','+B12','+B12','+B12','-B12','-B12','-B12')

temp_12 <- c('4','4','4','4','4','4','12','12','12','12','12','12')


sample_treat_12 <- data.frame(id_12,B12_12, temp_12)
colnames(sample_treat_12) <- c("id_12", "B12", "temp")

protquota_data <- cbind(protquota_data, sample_treat_12) |> 
  mutate(temp = factor(temp, levels = c("4", "12")),
         B12 = factor(B12, 
                      levels = c("+B12", 
                                 "-B12")),
         Treatment = factor(paste(temp, ",", B12)),
         fg_prot_cell = ug_cell_protein * 1e9)




# Protein Boxplot (Fig S3) ---------------------------------------------------------

protquota_bp1 <- ggplot(
  protquota_data,
  aes(
    x = temp,
    y = fg_prot_cell,
    grouping = Treatment,
    fill = temp,
    alpha = B12
  )
) +
  geom_boxplot() +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_fill_manual(values = col_12, labels = c("4", "12")) +
  scale_alpha_manual(values = c(1, .25), labels = c("+", "-")) +
  labs(
    y = bquote('Total Protein (Femtograms cell' ^ -1 * ')') ,
    x = "Treatment",
    fill = 'Temperature (°C)',
    alpha = expression("B"[12])
  ) +
  guides(alpha = guide_legend(override.aes = list(
    fill = hcl(c(15, 195), 100, 0, alpha = c(0.25, 1)), colour = NA
  ))) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) 



# Overall Protein t-test  -------------------------------------------------

# T-test to compare temps (p = 0.009954)
x_prot4 <- protquota_data |>  
  filter(temp == "4") |> 
  pull(fg_prot_cell)

y_prot12 <- protquota_data |> 
  filter(temp == "12") |> 
  pull(fg_prot_cell)

t.test(x_prot4, y_prot12)


