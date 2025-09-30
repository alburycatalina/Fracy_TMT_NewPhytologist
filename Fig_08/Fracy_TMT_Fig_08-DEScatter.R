# This script creates figure 8

# load libraries
library(wesanderson)
library(ggrepel)
library(ggpubr)


# set directory location 
here::i_am("Fig_08/Fracy_TMT_Fig_08-DEScatter.R")

kegg_annotated <- read.csv(here("./Fig_03/Fig_03_output_tables/kegg_annotated_17052025.csv"))


# Plot for B12 ------------------------------------------------------------

# Create a table with annotated proteins
B12_de_annotated <- kegg_annotated |>  filter(DE_origin == "B12")


# Group by A-level annotations
B12_de_annotated_A <- B12_de_annotated |> 
  group_by(accession, A) |> 
  dplyr::summarise(logFC= first(logFC), logCPM = first(logCPM))

# Grab a list of unique accession numbers
b12_accessions <- unique(B12_de_annotated$accession)

# List of names used in plot
b12_names <- c(
  "MetE",
  "protoporphyrin IX Mg-chelatase subunit D",
  "Unknown Protein (ID: 246327)",
  "P-ATPase family transporter",
  "T-complex protein 1 subunit gamma",
  "P-loop containing nucleoside triphosphate hydrolase protein",
  "ribosomal protein S12",
  "Unknown Protein (ID: 208848)",
  "Cys_Met_Meta_PP-domain-containing protein",
  "photosystem I p700 chlorophyll A apoprotein B",
  "Unknown Protein (ID: 260397)",
  "casein kinase I delta" ,
  "putative glucose-6-phosphate dehydrogenase" ,
  "porphobilinogen synthase" ,
  "Unknown Protein (ID: 271162)",
  "argininosuccinate synthase" ,
  "acetate CoA ligase",
  "TPR-like protein",
  "TPR-like protein",
  "Chloroa_b-bind-domain-containing protein",
  "mitochondrial 2-oxoglutarate/malate carrier protein",
  "phosphofructokinase",
  "chlorophyll a/b-binding protein",
  "cell division protease ftsH 4" ,
  "3-phosphoshikimate 1-carboxyvinyltransferase"
)

# Join in dataframe
b12_labels_df <- as.data.frame(cbind(b12_accessions, b12_names)) |> 
  rename("accession" = "b12_accessions",
         "label" = "b12_names")


b12_de_annotated_labeled <-
  B12_de_annotated_A |> 
  left_join(b12_labels_df, by = "accession")

# add chld 
g1 = subset(b12_de_annotated_labeled, accession == "OEU10229.1")

g2 = subset(b12_de_annotated_labeled, accession == "OEU11144.1")


b12_scatter <- ggplot(b12_de_annotated_labeled, aes(x = logCPM, y = logFC, color =
                                                      A)) +
  geom_point(size = 5, alpha = .50) +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  xlab(expression('Normalized Abundance')) +
  ylab(expression('Log'[2] * ' Fold Change')) +
  scale_color_manual(values = wes_palette("Royal2"), name = "Pathway") +
  geom_text_repel(data = g1, 
                  label = "ChlD2", 
                  vjust = -4, 
                  show.legend = FALSE, 
                  size = 5, 
                  nudge_x = 1.5)  + 
  geom_text_repel(data = g2, 
                  label = "MetE", 
                  # vjust = 1, 
                  show.legend = FALSE, 
                  size = 5, 
                  nudge_x = 1.5) +
  ylim(-7, 7) +
  xlim(8, 13) +
  theme(plot.title = element_text(hjust=0.5),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  ) +
  ggtitle(expression("-B"[12]))


# Plot for temp ------------------------------------------------------------

# Grab temp de's 
temp_de_annotated <- kegg_annotated |> 
  filter(DE_origin == "temp")


# List of labels for figure
temp_label <- c("Unknown Protein (ID: 271832)", "Unknown Protein (ID: 235337)")

# Temp scatter plot 
temp_scatter <- ggplot(temp_de_annotated, 
                       aes(x = logCPM, 
                           y = logFC, 
                           color = A)) +
  geom_point(size = 5, alpha = .50) +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  xlab(expression('Normalized Abundance')) +
  ylab(expression('Log'[2] * ' Fold Change')) +
  scale_color_manual(values = wes_palette("Royal2"), name = "Pathway") +
  geom_text_repel(
    aes(label = ifelse(logFC > 2.5 & logCPM > 8,
                       temp_label,
                       '')),
    
    hjust = 0,
    vjust = 0,
    show.legend =  FALSE,
    nudge_y = 1.5,
    size = 5
  ) +
  ylim(-7, 7) +
  ggtitle(expression("+12° C")) +
  theme(
    plot.title = element_text(hjust = 0.5),
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







# Plot for int DE proteins ------------------------------------------------
int_de_annotated <- kegg_annotated |>  filter(DE_origin == "int")


# Create table with shorter alt names for important proteins (logFC >/< 2.5; logCPM > 12)
accession <- c("OEU18748.1", "OEU11144.1", "OEU08040.1", "OEU11214.1", "OEU19288.1", "OEU15214.1", "OEU08856.1", "OEU22459.1", "OEU13453.1", "OEU15467.1", "OEU16276.1", "OEU11005.1")

names <- c("Unknown Protein (ID: 182880)", "MetE", "Unknown Protein (ID: 271832)", "CBA1", "Unknown Protein (ID: 235337)", "Nucleoside Triphosphate Hydrolase", "EngA", "NADH oxireductase", "Dihydrolipoamide acetyltransferase", "Triosephosphate isomerase", "Asparagine Synthase", "Band7 Domain Containing Protein" )

labels_df <- as.data.frame(cbind(accession, names))

int_de_annotated_labeled <- int_de_annotated |>  
  left_join(labels_df, by = "accession")

# Remove second MetE enrty to label
int_de_annotated_labeled <- int_de_annotated_labeled[-c(3),]

# Add ChlD for discussion purposes
g1 <- subset(int_de_annotated_labeled, accession == "OEU10229.1")

# Plotty plotty 
int_scatter <- ggplot(int_de_annotated_labeled, aes(x = logCPM, y = logFC, color =
                                                      A)) +
  geom_point(size = 5, alpha = .50) +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  xlab(expression('Normalized Abundance')) +
  ylab(expression('Log'[2] * ' Fold Change')) +
  scale_color_manual(values = wes_palette("Royal2"), name = "Pathway") +
  geom_text_repel(
    aes(label = ifelse(
      logFC < -2.5 & logCPM > 8, as.character(names), ''
    )),
    #hjust=0,
    #   vjust=0,
    show.legend =  FALSE,
    nudge_y = -1,
    nudge_x = 1,
    size = 5
  ) +
  geom_text_repel(
    aes(label = ifelse(
      logFC > 2.5 & logCPM > 8, as.character(names), ''
    )),
    #        hjust=0,
    #vjust=0,
    show.legend =  FALSE,
    nudge_y = 1,
    nudge_x = 1.75,
    size = 5
  ) +
  geom_text_repel(
    data = g1,
    label = "ChlD2",
    vjust = -1.5,
    hjust = -.75,
    show.legend = FALSE,
    size = 5,
    nudge_x = 1.5
  ) +
  xlim(7, 20) +
  ylim(-7, 7) +
  ggtitle(expression("-B"[12] * " and +12° C")) +
  theme(
    plot.title = element_text(hjust = 0.5),
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


plot <- ggarrange(b12_scatter, 
                  temp_scatter, 
                  int_scatter, 
                  ncol = 3, 
                  nrow = 1, 
                  common.legend = TRUE, 
                  legend = "bottom")

ggsave(plot, 
       filename = here("Fig_08/Fig_08_DEScatter.png"),  
       bg = "transparent", 
       width = 16, 
       height= 7, 
       units = "in")

