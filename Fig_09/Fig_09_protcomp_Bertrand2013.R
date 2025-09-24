
# This script makes figures that compare the MetE/MetH protein quotas to Bertrand et al. 2013


# Load required libraries 
library(tidyverse)
library(here)
library(MetBrewer)

# Load in and clean quota data 
fracy_metemeth_quotas <- read.csv(here("Fig_02/Fig_02_protein_quota_data/targeted_dataraw_treatsum.csv")) |> 
  
  # get only MetE and MetH
  filter(Temperature == "4",
         Protein != c("ACT1", "RBCL")) |> 
  # convert pmol to fmol 
  mutate(Treatment = "F. cylindrus (This study)",
         Organism = "F. cylindrus",
         Source = "This study",
         meanfmolAnalyte_ugProtein = treatmean_pmolAnalyte_ugProtein * 1000,
         sdfmolAnalyte_ugProtein = treatsd_pmolAnalyte_ugProtein * 1000) |> 
  select(
    Treatment,
    Organism,
    Source,
    Protein,
    B12,
    meanfmolAnalyte_ugProtein,
    sdfmolAnalyte_ugProtein
  )

# Load in data from Bertrand et al 2013 (P. trich and environmental)
Bertrand2013_metemeth_quotas <- read.csv(here("Fig_09/Bertrand2013_protquotas.csv"))

# Bind two df's 
prot_comp_df <- rbind(fracy_metemeth_quotas, Bertrand2013_metemeth_quotas) |> 
  mutate(Protein = factor(Protein,
         levels = c("MetH", "MetE")))


write.csv(prot_comp_df, file = here("Fig_09/prot_quota_comp.csv"))

# color palette for plot
b12_palette <- met.brewer("Hokusai1", 5)
b12_palette <- c(b12_palette[5], b12_palette[4], "light grey")

# Labels for plotting
fcyl_thisstudylabel <- ~ atop(paste(italic("F. cylindrus")), paste("(This study; 4° C)"))
pt_Bertrandlabel <- ~ atop(paste(italic("P.tricornutum")), paste("(Bertrand et al. 2013)"))
env_Bertrandlabel <- ~ atop(paste("McMurdo Sound"), paste("(Bertrand et al. 2013)"))


# Plot

prot_comp_plot <- ggplot(
  data = prot_comp_df,
  aes(
    y = meanfmolAnalyte_ugProtein,
    x = Treatment,
    group = Protein,
    color = Protein,
    shape = B12
  )
) +
  geom_point(size = 8, position = position_dodge(width = 0.5)) +
  ylab(expression(paste("Femtomoles µg total protein" ^ "-1"))) +
  guides(shape = guide_legend(
    title = expression("B"[12] ~ "Treatment"),
    override.aes = list(size = 4)
  )) +
  theme_classic() +
  xlab(NULL) +
  theme(text = element_text(size = 30),
        axis.text.x = element_text(
          size = 20,
          angle = 45,
          hjust = 1
        )) +
  scale_x_discrete(labels = c(fcyl_thisstudylabel, env_Bertrandlabel, pt_Bertrandlabel)) +
  scale_color_manual(values = b12_palette) +
  geom_errorbar(
    data = prot_comp_df,
    mapping = aes(ymin = meanfmolAnalyte_ugProtein - sdfmolAnalyte_ugProtein,
                  ymax = meanfmolAnalyte_ugProtein + sdfmolAnalyte_ugProtein),
    position = position_dodge(width = 0.5),
    width = .2
  ) +
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
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )

ggsave(prot_comp_plot, 
       file = here("Fig_09/Fig_09_protcomp_Bertrand2013.pdf"),  
       bg = "transparent",
       width = 14, 
       height = 9)

