

# Setup -------------------------------------------------------------------


#Load Packages
library(edgeR) # if installing edgeR, needs to be done through bioconductor. Instructions here: https://bioconductor.org/packages/release/bioc/html/edgeR.html
library(statmod)
library(tidyverse)

here::i_am("Fracy_TMT_DifferentialExpression/DE_20102021.R")


# Data cleaning -----------------------------------------------------------

# Load in file with protein quants and description
prot_data <- read.csv('Fracy_TMT_Normalization/PD1_Norm_11072021.csv')


# Select only the columns needed
DE_data <- prot_data |> select(
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

# Remove descriptions to include only counts
prot_counts <- DE_data[,2:length(DE_data)]

# Make vector with treatment groups
DE_data_groups <- factor(
  c(
    'B12_4',
    'B12_4',
    'B12_4',
    'noB12_4',
    'noB12_4',
    'B12_12',
    'B12_12',
    'noB12_12',
    'noB12_12'
  )
)

# Create DGEList object with protein counts, groups, and accession #
DGE_list <- DGEList(counts = prot_counts, 
                        group = DE_data_groups,
                        genes = DE_data$accession)

# Create design matrix 
design.mat <- model.matrix(~ 0 + DE_data_groups)

# Get tagwise dispersion of tags
DGE_disp <- estimateDisp(DGE_list, design.mat, robust = TRUE)


# Fit to GLM
fit <- glmQLFit(DGE_disp, design.mat, robust = TRUE)

# Make pairwise comparisons
# B12_4 vs B12_12
contrast_12_B12 <- makeContrasts(DE_data_groupsB12_12  - DE_data_groupsB12_4 , levels = design.mat)
# B12_4 vs noB12_4
contrast_4_noB12 <- makeContrasts(DE_data_groupsnoB12_4 - DE_data_groupsB12_4, levels = design.mat)
# B12_4 vs noB12_12
contrast_12_noB12 <- makeContrasts(  DE_data_groupsnoB12_12 - DE_data_groupsB12_4, levels = design.mat)

 # Fit tests

# B12_4 vs B12_12
qlf_12_B12 <- glmQLFTest(fit, contrast = contrast_12_B12)
# B12_4 vs noB12_4
qlf_4_noB12 <- glmQLFTest(fit, contrast = contrast_4_noB12)
# B12_4 vs noB12_12
qlf_12_noB12 <- glmQLFTest(fit, contrast = contrast_12_noB12)

# Export list of DE proteins per comparison

# B12_4 vs B12_12
hits_12_B12 <- topTags(qlf_12_B12, n = 500)
# B12_4 vs noB12_4
hits_4_noB12 <- topTags(qlf_4_noB12, n = 500)
# B12_4 vs noB12_12
hits_12_noB12 <- topTags(qlf_12_noB12, n = 500)

# Extract dataframes
hits_12_B12 <- hits_12_B12[[1]][,]
hits_4_noB12 <- hits_4_noB12[[1]][,]
hits_12_noB12 <- hits_12_noB12[[1]][,]

# Change genes column to say accession
hits_12_B12 <- hits_12_B12 |> 
  rename(accession = genes)
hits_4_noB12 <- hits_4_noB12 |> 
  rename(accession = genes)
hits_12_noB12 <- hits_12_noB12 |> 
  rename(accession = genes)

# Add back descriptions
hits_12_B12 <- merge(hits_12_B12, prot_data, by="accession")
hits_4_noB12 <- merge(hits_4_noB12, prot_data, by="accession")
hits_12_noB12 <- merge(hits_12_noB12, prot_data, by="accession")

# Trim dataset to only DE info
hits_12_B12 <- hits_12_B12 |> select(accession, Description, logFC, logCPM, F, PValue, FDR)
hits_4_noB12 <- hits_4_noB12 |> select(accession, Description, logFC, logCPM, F, PValue, FDR)
hits_12_noB12 <- hits_12_noB12 |> select(accession, Description, logFC, logCPM, F, PValue, FDR)

# Order by PValue 
hits_12_B12 <- hits_12_B12 |> arrange(PValue)
hits_4_noB12 <- hits_4_noB12 |> arrange(PValue)
hits_12_noB12 <- hits_12_noB12 |> arrange(PValue)


# Save hits lists to csv
write.csv(hits_12_B12, file = here('Fracy_TMT_DifferentialExpression/hits_12_B12_20102021.csv'))
write.csv(hits_4_noB12, file = here('Fracy_TMT_DifferentialExpression/hits_4_noB12_20102021.csv'))
write.csv(hits_12_noB12 , file = here('Fracy_TMT_DifferentialExpression/hits_12_noB12_20102021.csv'))

# Filter to only significant DE's (p< 0.05)
# Find out how many significant DE's
nDE_12_B12 <- nrow(hits_12_B12 |> filter(PValue < 0.05))
nDE_4_noB12 <- nrow(hits_4_noB12 |> filter(PValue < 0.05))
nDE_12_noB12 <- nrow(hits_12_noB12 |> filter(PValue < 0.05))

nDE <- c( nDE_12_B12, nDE_4_noB12, nDE_12_noB12)

