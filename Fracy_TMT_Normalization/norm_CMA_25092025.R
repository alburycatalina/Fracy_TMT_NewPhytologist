# Functions for TMT IRS Normalization with a variable number of MS experiments

# load required libraries
library(edgeR)
library(ggplot2)
library(dplyr)
library(magrittr)
library(testthat)
library(pastecs)
library(matrixStats)
library(devtools)
library(here)

# set directory location
here::i_am("norm_CMA_25092025.R")

#load and subset data
PD1 <- read.delim(here('190805_865_097_TMTfrag_plusPlastid_Proteins.txt'))

# Load sample names
PD1_names <- read.delim('new_col_names_11072021.txt') 
  


# re-order columns in dataset so that accession and quant are front and centre
PD1_sub <- PD1[c(4,19:38,2,5:18)]

# Use G-sub to replace column headers with info from col header file
names(PD1_sub) <- gsub(pattern = "Abundances.Grouped.",
                       replacement = "",
                       x = names(PD1_sub),
                       fixed = TRUE)


# Save new df with names added
PD1_sub2 <- PD1_sub


# Replace column headers
colnames(PD1_sub2) <- colnames(PD1_names)


# Re-order columns to give: accession, TMTsetA, TMTsetB, others. empty channels removed

PD1_sub3 <- PD1_sub2 [c(1,2,4,6,8,10,12,14,20,3,5,7,9,11,13,15,17,19,21,22:36)]


# Replace zeros in data with NAs

PD1_sub4A <- PD1_sub3[c(2:19)]
PD1_sub4B <- PD1_sub3[c(1,20:34)]

PD1_sub4A[PD1_sub4A == 0] <- NA

# Add back in columns from PD1_sub3 df
PD1_sub4 <- cbind(PD1_sub4A, PD1_sub4B)


PD1_sub4 <- PD1_sub4[c(19,1:18,20:34)]

colnames(PD1_sub4)

# basic stats for raw table
write.csv(PD1_sub4, file = "PD1_raw.csv")


options(scipen=100)
options(digits=2)
stats_PD1_sub4 <- stat.desc(PD1_sub4)

write.csv(stats_PD1_sub4, file = "stats_PD1_raw.csv")

# Remove low count experiment, 10 deg treatment, and outlier identified in cluster analysis noB12_12_1_B, B12_12_2_B
PD1_sub4_filt <- subset(PD1_sub4, select = -c(B12_10_1_A, noB12_10_1_A, B12_10_2_A, noB12_10_2_B ))

# Remove entries from crap database
PD1_sub4_filt  <-PD1_sub4_filt[!grepl("cRAP", PD1_sub4_filt$accession),]

# to impute rows with greater than n NAs. 
# Will choose 7 for FC (>1/2)
delete.na <- function(DF, n=0) {
  DF[rowSums(is.na(DF)) <= n,]
}


PD1_sub4_filt_7nas <- delete.na(PD1_sub4_filt, 7)


# Write data after removing rows with > 7 NA's
write.csv(PD1_sub4_filt_7nas, file = "PD1_sub4_filt_7nas.csv")



#imputing on NAs with row 1/2 min create fn

replace_na <- function(df){
  for(unique_row in 1:nrow(df)){
    # unique_row <- 1
    sub_df <- df[unique_row,]
    df[unique_row, is.na(sub_df)] <- sub_df$half_rowMin
  }
  return(df)
}  


# Imputation
# Get only protein values
PD1_sub5A_filt_7nas <- PD1_sub4_filt_7nas |> select(B12_4_1_A, noB12_4_1_A, B12_12_3_A, noB12_12_3_A, B12_12_1_B, noB12_12_2_B, B12_4_2_B, noB12_4_2_B, noB12_4_2_B, B12_4_3_B)


# Add a column that gives rowmin/2
PD1_sub5A_filt_7nas$rowMin <- apply(PD1_sub5A_filt_7nas, 1, FUN=min, na.rm = TRUE)
PD1_sub5A_filt_7nas$half_rowMin <- PD1_sub5A_filt_7nas$rowMin/2

# isolate accessions and other text from prev df
PD1_sub5B_filt_7nas <- PD1_sub4_filt_7nas |>  select(!c(B12_4_1_A, noB12_4_1_A, B12_12_3_A, noB12_12_3_A, B12_12_1_B, noB12_12_2_B, B12_4_2_B, noB12_4_2_B, noB12_4_2_B, B12_4_3_B))

# Replaced any NA's with half of row min
PD1_sub5A_filt_7imp <- replace_na(PD1_sub5A_filt_7nas)

# Add back in accessions and other text
PD1_sub5_filt_7imp <- cbind(PD1_sub5A_filt_7imp, PD1_sub5B_filt_7nas)



## check all rows are complete
PD1_sub6_filt_7imp <- PD1_sub5_filt_7imp[complete.cases(PD1_sub5_filt_7imp),]
write.csv(PD1_sub6_filt_7imp, file = "PD1_sub6_filt_7imp.csv")


# sample loading normalization
sl_normalization <- function(protein_df, tmt_exp_columns){
  
  # protein_df <- pd_light_all
  # tmt_exp_columns <- list(c(4:11), c(1:3))
  
  expect_is(tmt_exp_columns, 'list')
  expect_is(protein_df, 'data.frame')
  
  number_tmt_exps <- length(tmt_exp_columns)
  
  separated_prot_vals <- list()
  
  for(i in 1:number_tmt_exps){
    separated_prot_vals[[i]] <- protein_df[,tmt_exp_columns[[i]]]
  }
  
  # exp1_vals <- protein_df[,tmt1]
  # exp2_vals <- protein_df[,tmt2]
  
  normalization_factor_list <- list()
  
  for(i in 1:number_tmt_exps){
    normalization_factor_list[[i]] <- mean(colSums(separated_prot_vals[[i]])) / colSums(separated_prot_vals[[i]])
  }
  
  sl_normalized_prot_list <- list()
  
  for(i in 1:number_tmt_exps){
    sl_normalized_prot_list[[i]] <- sweep(separated_prot_vals[[i]], 2, normalization_factor_list[[i]], FUN = "*")
  }
  
  # col_blank <- rep('firebrick', length(box_labels))
  # col_blank[which(grepl(pattern = "_2$", x = box_labels))] <- 'darkblue'
  
  plot_df <- do.call('cbind', sl_normalized_prot_list)
  
  par(mar = c(7, 5, 3, 3))
  boxplot(log2(plot_df),
          # col = col_blank,
          xaxt = 'n',
          xlab = '',
          main = 'SL Normalization',
          ylab = 'Intensity')
  axis(1,
       labels = FALSE)
  text(x = seq_along(names(plot_df)),
       y = par("usr")[3] - 0.5,
       srt = 40,
       adj = 0.75,
       cex = 0.8,
       labels = names(plot_df),
       xpd = TRUE)
  
  return(sl_normalized_prot_list)
  
}

# Internal reference standard normalization
irs_normalization <- function(sl_normalized_list, tmt_common_channel_names){
  
  number_tmt_exps <- length(sl_normalized_list)
  
  irs_list <- list()
  irs_rowmeans_matrix <- matrix(nrow = length(rowMeans(sl_normalized_list[[1]])))
  
  for(i in 1:number_tmt_exps){
    # i <- 1
    sl_normalized_df <- sl_normalized_list[[i]]
    print( sl_normalized_df)
    print(tmt_common_channel_names[[i]])
    print(i)
    irs_list[[i]] <- sl_normalized_df[tmt_common_channel_names[[i]]]
    irs_rowmeans_matrix <- cbind(irs_rowmeans_matrix, rowMeans(irs_list[[i]]))
  }
  
  rowmean_vector <- rowMeans(irs_rowmeans_matrix, na.rm = TRUE)
  
  scaling_factor_list <- list()
  irs_sl_normalized_list <- list()
  
  for(i in 1:number_tmt_exps){
    scaling_factor_list[[i]] <- rowmean_vector / rowMeans(irs_list[[i]])
    sl_normalized_df <- sl_normalized_list[[i]]
    irs_sl_normalized_list[[i]] <- sl_normalized_df * scaling_factor_list[[i]]
  }
  
  irs_sl_df <- matrix(nrow = nrow(irs_sl_normalized_list[[1]])) |> as.data.frame()
  
  for(i in 1:number_tmt_exps){
    irs_sl_df <- cbind(irs_sl_df, irs_sl_normalized_list[[i]])
  }
  
  return_df <- irs_sl_df[,-1]
  # 
  # col_blank <- rep('firebrick', length(box_labels))
  # col_blank[which(grepl(pattern = "_2$", x = box_labels))] <- 'darkblue'
  # 
  par(mar = c(7, 5, 3, 3))
  boxplot(log2(return_df),
          # col = col_blank,
          xaxt = 'n',
          xlab = '',
          main = 'SL, IRS Normalization',
          ylab = 'Intensity')
  axis(1,
       labels = FALSE)
  text(x = seq_along(names(return_df)),
       y = par("usr")[3] - 0.5,
       srt = 40,
       adj = 0.75,
       cex = 0.8,
       labels = names(return_df),
       xpd = TRUE)
  
  return(return_df)
  
}


# edgeR tmm normalization
tmm_normalization <- function(irs_normalized_df){
  
  irs_tmm <- edgeR::calcNormFactors(irs_normalized_df)
  data_irs_tmm <- sweep(irs_normalized_df, 2, irs_tmm, FUN = "/")
  
  
  boxplot(log2(data_irs_tmm), 
          # col = c(rep('firebrick', length(tmt1)),
          # rep('darkblue', length(tmt2))),
          # col = col_blank,
          xaxt = 'n', 
          xlab = '',
          ylab = 'Intensity',
          main = 'SL, IRS, TMM Normalization')
  axis(1, 
       labels = FALSE)
  text(x = seq_along(names(data_irs_tmm)), 
       y = par("usr")[3] - 0.5, 
       srt = 45, 
       adj = 1,
       labels = names(data_irs_tmm), 
       xpd = TRUE)
  
  return(data_irs_tmm)
}


sl_irs_tmm_normalization <- function(protein_df, tmt_exp_columns, tmt_common_channel_names){
  
  # check that the inputs are of the correct type
  expect_is(protein_df, 'data.frame')
  expect_is(tmt_exp_columns, 'list')
  expect_is(tmt_common_channel_names, 'list')
  
  # checking input files
  is.contained <- function(x, y) {
    z <- x[x %in% setdiff(x, y)]
    length(z) == length(x) - length(y)
  }
  
  # checking that the common channel names are within the protein df
  testing_tmt_names <- is.contained(x = names(protein_df), y = unlist(tmt_common_channel_names))
  
  if(!testing_tmt_names){
    stop('It seems like the tmt_common_channel_names contains a name that is not in your protein_df; double check that!')
  }
  
  complete_protein_df <- protein_df[complete.cases(protein_df),]
  
  if(nrow(complete_protein_df) != nrow(protein_df)){
    warning('Your protein_df has some missing values. This normalization can only consider proteins observed across all TMT channels. The output of this function subsets only the normalized proteins observed across all channels!')
    protein_df <- complete_protein_df
  }
  
  for(i in 1:length(tmt_exp_columns)){
    if(i == 1) print('These are the tmt experiment column labels youve designated, just to be sure:')
    print('---------------------------------')
    print(paste0('TMT Experiment ', i))
    names(protein_df)[tmt_exp_columns[[i]]] |> print()
    print('---------------------------------')
  }
  
  par(mfrow = c(2, 2))
  
  boxplot(log2(protein_df), 
          # col = c(rep('firebrick', length(tmt1)),
          # rep('darkblue', length(tmt2))),
          # col = col_blank,
          xaxt = 'n', 
          xlab = '',
          ylab = 'Intensity',
          main = 'No Normalization')
  axis(1, 
       labels = FALSE)
  text(x = seq_along(names(protein_df)), 
       y = par("usr")[3] - 0.5, 
       srt = 45, 
       adj = 1,
       labels = names(protein_df), 
       xpd = TRUE)
  sl_norm_file <- sl_normalization(protein_df = protein_df, tmt_exp_columns = tmt_exp_columns)
  irs_norm_file <- irs_normalization(sl_normalized_list = sl_norm_file, tmt_common_channel_names = tmt_common_channel_names)
  tmm_norm_file <- tmm_normalization(irs_normalized_df = irs_norm_file)
  return(tmm_norm_file)
}

# Get a df with only sample and pool values
PD1_sub6_filt_7imp_just_val <- PD1_sub6_filt_7imp |> select(
  B12_4_1_A,
  noB12_4_1_A,
  B12_12_3_A,
  noB12_12_3_A,
  B12_12_1_B,
  noB12_12_2_B,
  B12_4_2_B,
  noB12_4_2_B,
  noB12_4_2_B,
  B12_4_3_B,
  pool2_A,
  pool1_B,
  pool2_B
)

# Pools are in 
# Running the normalization:
PD1_Norm <- sl_irs_tmm_normalization(
  protein_df = PD1_sub6_filt_7imp_just_val,
  # list experiment including pools by experiment (A then B)
  tmt_exp_columns = list(c(1:4, 10), c(5:9, 11:12)),
  #give pool names
  tmt_common_channel_names = list(c("pool2_A"), c("pool1_B", "pool2_B"))
)

pdf(file = here("../Fig_S4/Fig_S4.pdf"), width = 12, height = 6)

dev.off()

# Add back in accessions etc. 
PD1_Norm_final <-  cbind(
  PD1_sub6_filt_7imp |> select(
    accession,
    Description,
    Protein.FDR.Confidence.Combined,
    Exp.q.value.Combined,
    Sum.PEP.Score,
    Coverage.in.Percent,
    Number.of.Peptides,
    Number.of.PSMs,
    Number.of.Unique.Peptides
  ) ,
  PD1_Norm
) |> write.csv(file = "PD1_Norm_11072021.csv")

