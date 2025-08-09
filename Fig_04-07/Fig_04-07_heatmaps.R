# This script completes a hierarchical clustering analysis and plots heatmaps from differential expression on proteins from Fracy TMT multiplex experiment 



# Setup -------------------------------------------------------------------


# Load required packages
library(dplyr)
library(gplots)
library(stringr)
library(pheatmap)
library(sjmisc)
library(here)

here::i_am("Fig_04-07/Fig_04-07_heatmaps.R")


# Clean up data -----------------------------------------------------------


# Load entire normalized dataset; change class of description column to matrix
norm_data <- read.csv(here('Fracy_TMT_Normalization/PD1_Norm_11072021.csv')) |>
  mutate(Description = make.unique(as.character(Description))) |>
  select( # include only certain columns with accession, description, and protein expression data (excludes 10deg C treatments and the low count noB12, 12 treatment)
    accession,
    B12_4_1_A,
    B12_4_2_B,
    B12_4_3_B,
    noB12_4_1_A,
    noB12_4_2_B,
    B12_12_1_B,
    B12_12_3_A,
    noB12_12_2_B,
    noB12_12_3_A
  )


# Load 3 toptags B12 hits from DE analysis pairwise comparisons
hits_4_noB12  <- read.csv(here('Fracy_TMT_DifferentialExpression/hits_4_noB12_20102021.csv'))
hits_12_B12 <- read.csv(here('Fracy_TMT_DifferentialExpression/hits_12_B12_20102021.csv'))
hits_12_noB12  <- read.csv(here('Fracy_TMT_DifferentialExpression/hits_12_noB12_20102021.csv'))



# B12 Heatmap -------------------------------------------------------------

# Filter out significant B12 hits
hits_4_noB12_sigonly  <- hits_4_noB12 |> filter(PValue < .01) |> 
  # change accession cols to character type
  mutate(accesssion = as.character(accession)) 

# Match hits to normalized data by accession #'s
hits_sig_joined_B12 <- inner_join(hits_4_noB12_sigonly,  norm_data,
                                  by=c("accession" = "accession"))

# Remove unwanted col's (anything but count data)
heatmap_hits_B12 <- as.matrix(
  dplyr::select(
    hits_sig_joined_B12,
    B12_4_1_A,
    B12_4_2_B,
    B12_4_3_B,
    B12_12_1_B,
    B12_12_3_A,
    noB12_4_1_A,
    noB12_4_2_B,
    noB12_12_2_B,
    noB12_12_3_A
  )
)


# Make row names from accession column in hits_sig_joined
rownames(heatmap_hits_B12) = make.names(hits_sig_joined_B12$accession, unique=TRUE)



# Scale heatmap hits to remove long tail issues (https://stackoverflow.com/questions/21983162/how-to-expand-the-dendogram-in-heatmap-2)
heatmap_hits_scaled_B12 <- t(scale(t(heatmap_hits_B12))) 

# set custom distance and clustering functions
hclustfunc <- function(x) hclust(x, method="average")
distfunc <- function(x) dist(x,method="maximum")


# Create df for annotations
annotation_row_B12 <- data.frame(hits_sig_joined_B12$Description)
rownames(annotation_row_B12) <- rownames(heatmap_hits_scaled_B12)

# Save hits_sig_joined
write.csv(hits_sig_joined_B12, file = here("Fig_04-07/Fig_04-07_OutputTables/de_heatmap_hits_B12.csv"))

# Pairwise correlation between columns (treatments)
cols.cor <- cor(heatmap_hits_scaled_B12, use = "pairwise.complete.obs", method = "pearson")
# Pairwise correlation between rows (proteins)
rows.cor <- cor(t(heatmap_hits_scaled_B12), use = "pairwise.complete.obs", method = "pearson")


annotations <- read.csv(here("Fig_04-07/Fig_04-07_AnnotationTables/de_heatmap_hits_B12_annotated.csv"))
annotation_row_B12 <- data.frame(annotations$accession)
rownames(annotation_row_B12) <- rownames(heatmap_hits_scaled_B12)
colnames(annotation_row_B12)[1] <- "Function"

# Change order
col.order <- c("B12_4_1_A", "B12_4_2_B", "B12_4_3_B", "noB12_4_1_A", "noB12_4_2_B", "B12_12_1_B", "B12_12_3_A", "noB12_12_2_B", "noB12_12_3_A")

heatmap_hits_scaled_B12 <- heatmap_hits_scaled_B12[,col.order]

# Plot the heatmap
b12_heatmap <- pheatmap(heatmap_hits_scaled_B12,
         scale = "row", 
         cluster_cols = FALSE, 
         treeheight_row = 150,
         clustering_distance_rows = as.dist(1 - rows.cor),
         legend = TRUE, 
         labels_row = annotations$heatmap_rowname, 
         labels_col = c("4, +B12", 
                        "4, +B12", 
                        "4, +B12", 
                        "4, -B12",
                        "4, -B12",
                        "12, +B12",
                        "12, +B12", 
                        "12, -B12", 
                        "12, -B12"), 

         fontsize_row = 13, 
         fontsize_col = 17,
         border_color = "black",
         gaps_col = c(3,5,7))

ggsave(b12_heatmap,
       file = here("Fig_04-07/Fig_04-07_HeatmapFigs/Fig_05_B12heatmap.pdf"), 
       width = 12,
       height = 8)


# Heatmap for temp comp ---------------------------------------------------


# Filter out significant B12 hits
hits_12_B12_sigonly  <- hits_12_B12 |> filter(PValue < .01)

# Change accession cols to character type
hits_12_B12_sigonly$accession <- as.character(hits_12_B12_sigonly$accession)

# Match hits to normalized data by accession #'s
hits_sig_joined_temp <- inner_join(hits_12_B12_sigonly,  norm_data,
                                  by=c("accession" = "accession"))

# Remove unwanted col's (anything but count data)
heatmap_hits_temp <- as.matrix(select(hits_sig_joined_temp, B12_4_1_A, B12_4_2_B, B12_4_3_B, B12_12_1_B, B12_12_3_A, noB12_4_1_A, noB12_4_2_B, noB12_12_2_B, noB12_12_3_A))


# Make row names from accession column in hits_sig_joined
rownames(heatmap_hits_temp) = make.names(hits_sig_joined_temp$accession, unique=TRUE)

# Make heatmap Row names
# Add MetE, CBA1, thiC, nucleoside triphosphate hydrolase protein entries
hits_sig_joined_temp$Description[hits_sig_joined_temp$accession == "OEU11144.1"] <- "MetE"
hits_sig_joined_temp$Description[hits_sig_joined_temp$accession == "OEU16390.1"] <- "ThiC"
hits_sig_joined_temp$Description[hits_sig_joined_temp$accession == "OEU11214.1"] <- "CBA1"
hits_sig_joined_temp$Description[hits_sig_joined_temp$accession == "OEU15214.1"] <- "Nucleoside Triphosphate Hydrolase"
hits_sig_joined_temp$Description[hits_sig_joined_temp$accession == "OEU18445.1"] <- "Nucleoside Triphosphate Hydrolase"
hits_sig_joined_temp$Description[hits_sig_joined_temp$accession == "OEU23066.1"] <- "Nucleoside Triphosphate Hydrolase"


colnames(hits_sig_joined_temp)[3] <- "Description"

# Remove f. cylindrus from end of protein names
hits_sig_joined_temp$Description <- str_remove(hits_sig_joined_temp$Description, fixed(pattern = "[Fragilariopsis cylindrus CCMP1102]"))


for (i in 1:nrow(hits_sig_joined_temp)){
  if (str_contains(hits_sig_joined_temp$Description[i], "hypothetical")){
    # Add column with protein id's for hypothetical proteins
    hits_sig_joined_temp$proteinId[i] <- gsub(".*_", "", hits_sig_joined_temp$Description[i]) 
    hits_sig_joined_temp$proteinId[i] <- gsub(" ", "", hits_sig_joined_temp$proteinId[i], fixed = TRUE)
    
    hits_sig_joined_temp$heatmap_rowname[i] <- paste0("Uncharacterized Protein (ID:", hits_sig_joined_temp$proteinId[i],")")
  }
  else{
    hits_sig_joined_temp$heatmap_rowname[i] <- hits_sig_joined_temp$Description[i]
  }
}


# Scale heatmap hits to remove long tail issues (https://stackoverflow.com/questions/21983162/how-to-expand-the-dendogram-in-heatmap-2)

heatmap_hits_scaled_temp <- t(scale(t(heatmap_hits_temp))) 

# set custom distance and clustering functions
hclustfunc <- function(x) hclust(x, method="average")
distfunc <- function(x) dist(x,method="maximum")


# Create df for annotations
annotation_row_temp <- data.frame(hits_sig_joined_temp$Description)
rownames(annotation_row_temp) <- rownames(heatmap_hits_scaled_temp)

# Save hits_sig_joined
write.csv(hits_sig_joined_temp, file = here("Fig_04-07/Fig_04-07_OutputTables/de_heatmap_hits_temp.csv"))

# Pairwise correlation between columns (treatments)
cols.cor <- cor(heatmap_hits_scaled_temp, use = "pairwise.complete.obs", method = "pearson")
# Pairwise correlation between rows (proteins)
rows.cor <- cor(t(heatmap_hits_scaled_temp), use = "pairwise.complete.obs", method = "pearson")


# annotations_temp <- read.csv("de_heatmap_hits_temp_annotated.csv")
# annotation_row_temp <- data.frame(annotations_temp$accession)
# rownames(annotation_row_temp) <- rownames(heatmap_hits_scaled_temp)
# colnames(annotation_row_temp)[1] <- "Function"

# Change order
col.order <- c("B12_4_1_A", "B12_4_2_B", "B12_4_3_B", "noB12_4_1_A", "noB12_4_2_B", "B12_12_1_B", "B12_12_3_A", "noB12_12_2_B", "noB12_12_3_A")

heatmap_hits_scaled_temp <- heatmap_hits_scaled_temp[,col.order]

# Plot the heatmap
temp_heatmap <- pheatmap(heatmap_hits_scaled_temp,
         scale = "row", 
         cluster_cols = FALSE, 
         treeheight_row = 150,
         clustering_distance_rows = as.dist(1 - rows.cor),
         legend = TRUE, 
         labels_row = hits_sig_joined_temp$heatmap_rowname, 
         labels_col = c("4, +B12", 
                        "4, +B12", 
                        "4, +B12", 
                        "4, -B12",
                        "4, -B12",
                        "12, +B12",
                        "12, +B12", 
                        "12, -B12", 
                        "12, -B12"), 
         border_color = "black",
         fontsize_row = 11, 
         fontsize_col = 15,
         gaps_col = c(3,5,7)
)

ggsave(temp_heatmap,
       file = here("Fig_04-07/Fig_04-07_HeatmapFigs/Fig_06_tempheatmap.pdf"), 
       width = 12,
       height = 8)




# Heatmap for all interaction comp ----------------------------------------


# Filter out significant B12 hits
hits_12_noB12_sigonly  <- hits_12_noB12 |> filter(PValue < .01)

# Change accession cols to character type
hits_12_noB12_sigonly$accession <- as.character(hits_12_noB12_sigonly$accession)

# Match hits to normalized data by accession #'s
hits_sig_joined_int <- inner_join(hits_12_noB12_sigonly,  norm_data,
                                   by=c("accession" = "accession"))

# Remove unwanted col's (anything but count data)
heatmap_hits_int <- as.matrix(select(hits_sig_joined_int, B12_4_1_A, B12_4_2_B, B12_4_3_B, B12_12_1_B, B12_12_3_A, noB12_4_1_A, noB12_4_2_B, noB12_12_2_B, noB12_12_3_A))


# Make row names from accession column in hits_sig_joined
rownames(heatmap_hits_int) = make.names(hits_sig_joined_int$accession, unique=TRUE)



# Scale heatmap hits to remove long tail issues (https://stackoverflow.com/questions/21983162/how-to-expand-the-dendogram-in-heatmap-2)

heatmap_hits_scaled_int <- t(scale(t(heatmap_hits_int))) 

# set custom distance and clustering functions
hclustfunc <- function(x) hclust(x, method="average")
distfunc <- function(x) dist(x,method="maximum")


# Create df for annotations
annotation_row_int <- data.frame(hits_sig_joined_int$Description)
rownames(annotation_row_int) <- rownames(heatmap_hits_scaled_int)

# Make heatmap Row names
# Add MetE, CBA1, thiC, nucleoside triphosphate hydrolase protein entries
hits_sig_joined_int$Description[hits_sig_joined_int$accession == "OEU11144.1"] <- "MetE"
hits_sig_joined_int$Description[hits_sig_joined_int$accession == "OEU16390.1"] <- "ThiC"
hits_sig_joined_int$Description[hits_sig_joined_int$accession == "OEU11214.1"] <- "CBA1"
hits_sig_joined_int$Description[hits_sig_joined_int$accession == "OEU15214.1"] <- "Nucleoside Triphosphate Hydrolase"
hits_sig_joined_int$Description[hits_sig_joined_int$accession == "OEU18445.1"] <- "Nucleoside Triphosphate Hydrolase"

colnames(hits_sig_joined_int)[3] <- "Description"

# Remove f. cylindrus from end of protein names
hits_sig_joined_int$Description <- str_remove(hits_sig_joined_int$Description, fixed(pattern = "[Fragilariopsis cylindrus CCMP1102]"))


for (i in 1:nrow(hits_sig_joined_int)){
  if (str_contains(hits_sig_joined_int$Description[i], "hypothetical")){
    # Add column with protein id's for hypothetical proteins
    hits_sig_joined_int$proteinId[i] <- gsub(".*_", "", hits_sig_joined_int$Description[i]) 
    hits_sig_joined_int$proteinId[i] <- gsub(" ", "", hits_sig_joined_int$proteinId[i], fixed = TRUE)

    hits_sig_joined_int$heatmap_rowname[i] <- paste0("Uncharacterized Protein (ID:", hits_sig_joined_int$proteinId[i],")")
  }
else{
  hits_sig_joined_int$heatmap_rowname[i] <- hits_sig_joined_int$Description[i]
}
}
# save df
write.csv(hits_sig_joined_int, file = here("Fig_04-07/Fig_04-07_OutputTables/de_heatmap_hits_int.csv"))

# Pairwise correlation between columns (treatments)
cols.cor <- cor(heatmap_hits_scaled_int, use = "pairwise.complete.obs", method = "pearson")
# Pairwise correlation between rows (proteins)
rows.cor <- cor(t(heatmap_hits_scaled_int), use = "pairwise.complete.obs", method = "pearson")


#annotations_int <- read.csv("de_heatmap_hits_temp_annotated.csv")
#annotation_row_int <- data.frame(hits_sig_joined_int$accession)
#rownames(annotation_row_int) <- rownames(heatmap_hits_scaled_int)
#colnames(annotation_row_int)[1] <- "Function"

# Change order
col.order <- c("B12_4_1_A", "B12_4_2_B", "B12_4_3_B", "noB12_4_1_A", "noB12_4_2_B", "B12_12_1_B", "B12_12_3_A", "noB12_12_2_B", "noB12_12_3_A")

heatmap_hits_scaled_int <- heatmap_hits_scaled_int[,col.order]

# Plot the heatmap
int_heatmap <- pheatmap(heatmap_hits_scaled_int,
         scale = "row", 
         cluster_cols = FALSE, 
         treeheight_row = 150,
         clustering_distance_rows = as.dist(1 - rows.cor),
         legend = TRUE, 
         labels_row = hits_sig_joined_int$heatmap_rowname, 
         labels_col = c("4, +B12", 
                        "4, +B12", 
                        "4, +B12", 
                        "4, -B12",
                        "4, -B12",
                        "12, +B12",
                        "12, +B12", 
                        "12, -B12", 
                        "12, -B12"), 
         border_color = "black",
         fontsize_row = 7, 
         fontsize_col = 15,
         gaps_col = c(3,5,7), 
         legend_labels = c("thing")
)

ggsave(int_heatmap,
       file = here("Fig_04-07/Fig_04-07_HeatmapFigs/Fig_07_intheatmap.pdf"), 
       width = 12,
       height = 8)




# Heatmap With All Pairwise Comps -----------------------------------------
# Filter toptags to only significantly DE proteins (p = 0.01 and log counts per million = 10.5)
hits_4_noB12_sig  <- hits_4_noB12 |> filter(PValue < .01 & logCPM > 10.5)
hits_12_B12_sig <- hits_12_B12 |> filter(PValue < .01 & logCPM > 10.5)
hits_12_noB12_sig  <- hits_12_noB12 |> filter(PValue < .01 & logCPM > 10.5)

# Change accession cols to character type
hits_4_noB12_sig$accession <- as.character(hits_4_noB12_sig$accession)
hits_12_B12_sig$accession <- as.character(hits_12_B12_sig$accession)
hits_12_noB12_sig$accession <- as.character(hits_12_noB12_sig$accession)

# Merge 3 data sets vertically so there are no duplicate protein names (Collapse by 'genes' col)
hits_merged <- rbind( hits_12_B12_sig, hits_4_noB12_sig, hits_12_noB12_sig)

# Match hits to normalized data by accession #'s
hits_sig_joined <- inner_join(hits_merged,  norm_data,
                              by=c("accession" = "accession"))

# Check for duplicate rows
hits_sig_joined$accession[duplicated(hits_sig_joined$accession)]


# There are some (same proteins are DE in multiple treatments) - remove
hits_sig_joined <- hits_sig_joined |> distinct(accession, .keep_all = TRUE)


# Remove unwanted col's (anything but count data)
heatmap_hits <- as.matrix(select(hits_sig_joined, B12_4_1_A, B12_4_2_B, B12_4_3_B, B12_12_1_B, B12_12_3_A, noB12_4_1_A, noB12_4_2_B, noB12_12_2_B, noB12_12_3_A))


# Make row names from accession column in hits_sig_joined
rownames(heatmap_hits) = make.names(hits_sig_joined$accession, unique=TRUE)

# Scale heatmap hits to remove long tail issues (https://stackoverflow.com/questions/21983162/how-to-expand-the-dendogram-in-heatmap-2)

heatmap_hits_scaled <- t(scale(t(heatmap_hits))) 

# set custom distance and clustering functions
hclustfunc <- function(x) hclust(x, method="average")
distfunc <- function(x) dist(x,method="maximum")

# Remove "[Fragilariopsis cylindrus CCMP1102]" from end of all descriptions
colnames(hits_sig_joined)[3] <- "Description"

hits_sig_joined$Description <- str_remove(hits_sig_joined$Description, fixed(pattern = "[Fragilariopsis cylindrus CCMP1102]"))


# Create df for annotations
annotation_row <- data.frame(hits_sig_joined$Description)
rownames(annotation_row) <- rownames(heatmap_hits_scaled)

# Save hits_sig_joined
write.csv(hits_sig_joined, file = here("Fig_04-07/Fig_04-07_OutputTables/de_heatmap_hits.csv"))

# Pairwise correlation between columns (treatments)
cols.cor <- cor(heatmap_hits_scaled, use = "pairwise.complete.obs", method = "pearson")
# Pairwise correlation between rows (proteins)
rows.cor <- cor(t(heatmap_hits_scaled), use = "pairwise.complete.obs", method = "pearson")

# Change MetE description
hits_sig_joined$Description[hits_sig_joined$accession == "OEU11144.1"] <- "MetE"

annotations <- read.csv("Fig_04-07/Fig_04-07_AnnotationTables/de_heatmap_hits_annotated_11072021.csv")
annotation_row <- data.frame(annotations$Function)
rownames(annotation_row) <- rownames(heatmap_hits_scaled)
colnames(annotation_row)[1] <- "Function"

# Make heatmap Row names
# Add MetE, CBA1, thiC, nucleoside triphosphate hydrolase protein entries
hits_sig_joined$Description[hits_sig_joined$accession == "OEU11144.1"] <- "MetE"
hits_sig_joined$Description[hits_sig_joined$accession == "OEU16390.1"] <- "ThiC"
hits_sig_joined$Description[hits_sig_joined$accession == "OEU11214.1"] <- "CBA1"
hits_sig_joined$Description[hits_sig_joined$accession == "OEU15214.1"] <- "Nucleoside Triphosphate Hydrolase"
hits_sig_joined$Description[hits_sig_joined$accession == "OEU18445.1"] <- "Nucleoside Triphosphate Hydrolase"

colnames(hits_sig_joined)[3] <- "Description"

hits_sig_joined$proteinId <- NA

for (i in 1:nrow(hits_sig_joined)){
  if (str_contains(hits_sig_joined$Description[i], "hypothetical")){
    # Add column with protein id's for hypothetical proteins
    hits_sig_joined$proteinId[i] <- gsub(".*_", "", hits_sig_joined$Description[i])
    hits_sig_joined$proteinId[i] <- gsub(" ", "", hits_sig_joined$proteinId[i], fixed = TRUE)

    hits_sig_joined$heatmap_rowname[i] <- paste0("Uncharacterized Protein (ID:", hits_sig_joined$proteinId[i],")")
  }
  else{
    hits_sig_joined$heatmap_rowname[i] <- hits_sig_joined$Description[i]
  }
}

# Change order
col.order <- c("B12_4_1_A", "B12_4_2_B", "B12_4_3_B", "noB12_4_1_A", "noB12_4_2_B", "B12_12_1_B", "B12_12_3_A", "noB12_12_2_B", "noB12_12_3_A")

heatmap_hits_scaled <- heatmap_hits_scaled[,col.order]



# Plot the heatmap
heatmap_all <- pheatmap(heatmap_hits_scaled,
         scale = "row", 
         cluster_cols = FALSE, 
         treeheight_row = 150,
         clustering_distance_rows = as.dist(1 - rows.cor),
         legend = TRUE, 
         labels_row = hits_sig_joined$heatmap_rowname, 
         labels_col = c("4, +B12", 
                        "4, +B12", 
                        "4, +B12", 
                        "4, -B12",
                        "4, -B12",
                        "12, +B12",
                        "12, +B12", 
                        "12, -B12", 
                        "12, -B12"), 
         border_color = "black",
         fontsize_row = 11, 
         fontsize_col = 15,
         gaps_col = c(3,5,7)
)

ggsave(heatmap_all,
       file = here("Fig_04-07/Fig_04-07_HeatmapFigs//Fig_04_heatmap.pdf"), 
       width = 12,
       height = 8)

