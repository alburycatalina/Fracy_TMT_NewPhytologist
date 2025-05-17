# Script for functional annotation of DE proteins from TMT experiment with KEGG values


# Setup :-) ---------------------------------------------------------------

# Load required packages into environment
library(tidyverse)
library(KEGGREST)
library(reshape2)
library(ggrepel)
library(wesanderson)
library(ggpubr)
library(kableExtra)
library(stringr)
library(data.table)
library(grid)
library(scales)

# Installing KEGGREST
# BiocManager::install("KEGGREST")


# Load DE Data ------------------------------------------------------------


# Load in data with hits of DE'd proteins 4,B12 vs 4,noB12  (B12 response)
B12_hits <- read.csv(here("Fig_3/Fig_3_raw_data/hits_4_noB12_20102021.csv")) |>
  mutate(Description = str_remove(Description, fixed(pattern = "[Fragilariopsis cylindrus CCMP1102]")),
         DE_origin = "B12")



# 4,B12 vs 12,B12 (temp response)
temp_hits <- read.csv(here("Fig_3/Fig_3_raw_data/hits_12_B12_20102021.csv")) |>
  mutate(Description = str_remove(Description, fixed(pattern = "[Fragilariopsis cylindrus CCMP1102]")),
         DE_origin = "temp")

# 4,B12 vs 12,noB12 (interaction repsponse between B12 and temp)
int_hits <- read.csv(here("Fig_3/Fig_3_raw_data/hits_12_noB12_20102021.csv")) |>
mutate(Description = str_remove(Description, fixed(pattern = "[Fragilariopsis cylindrus CCMP1102]")),
       DE_origin = "int")


# Merge hit lists into 1, only including significantly DE'd proteins (Pvalues < .05)
hits_list <- rbind(B12_hits, temp_hits, int_hits) |> 
  filter(PValue < .05) |>
  # Make a column that indicates DE direction (up or downregulated; positive or negative fold change)
  mutate(de_dir = fifelse(logFC < 0,
                          -1,
                          1))

# Add KEGG Annotation Information -----------------------------------------

# Match in protein ID's from NCBI (https://www.ncbi.nlm.nih.gov/genome/browse/#!/proteins/11246/283952%7CFragilariopsis%20cylindrus%20CCMP1102/)
prot_info <- read.csv(here("Fig_3/Fig_3_raw_data/prot_table.csv")) |>
  rename("accession" = "Protein.product")

# Protein ID's for chloroplast from (https://www.ncbi.nlm.nih.gov/genome/browse/#!/proteins/11246/748136%7CFragilariopsis%20cylindrus/chloroplast/)
# chl_prot_info <- read.csv("chl_prot_table.csv")

# Match accessions to protein ID's
prot_data_ids <- left_join(hits_list, prot_info, by = "accession") |>
  
  # Make new column with just protein ID's (6 digit code at end of locus.tag)
  mutate(proteinId = gsub(".*_", "", Locus.tag))



# # Match ID's to functional information from KEGG ------------------------

# Use KEGGREST package to access KEGG REST API 
# https://www.bioconductor.org/packages/release/bioc/vignettes/KEGGREST/inst/doc/KEGGREST-vignette.html

# Save fracy data to rds file for quick loading
# saveRDS(as.data.frame(keggList("fcy")), 
        # here("Fig_3/Fig_3_raw_data", "KEGGlist_frag.rds"))

# Load in RDS if available in cache
KEGGlist_frag <- readRDS(here("Fig_3/Fig_3_raw_data", "KEGGlist_frag.rds"))



# FIXME source for these matching KEGG numbers?
frag_kegg_nos_raw <- read.csv(here("Fig_3/Fig_3_raw_data/frag_kegg_no.csv"))

# Parse columns by separators
frag_kegg_nos <- colsplit(frag_kegg_nos_raw$Name," ",c("id","kegg_id"))

temp <- colsplit(frag_kegg_nos$kegg_id," ", c("K_no","description"))
database <- cbind(frag_kegg_nos, temp)

# Delete row if no KO number is available
frag_kegg <- database[!grepl("no", database$K_no),] |> 
  
  # Remove unneeded col
  select(-c("kegg_id"))


# Produce a dataframe matching protein id's and kegg numbers
frag_kegg_database <- frag_kegg |> mutate(colsplit(id, 
                                            pattern = "_",
                                            names = c("fracy","proteinId"))) |> 
  select(K_no, proteinId)


# Join protein data and keg no's
frag_kegg_final <- left_join(prot_data_ids, frag_kegg_database, by = "proteinId")



# Tables of DE Proteins ---------------------------------------------------


# Create tables for paper
# B12 table 
B12_table <- filter(frag_kegg_final, DE_origin == "B12", PValue < .05) |>  
  select(accession, 
         Description, 
         proteinId, 
         K_no, 
         logFC, 
         logCPM, 
         PValue)

# Write to SI folder
write.csv(B12_table, 
          file = here("Fig_S5-S7/CMA_NewPhytologist_B12_table_TMT_17052025.csv"))


# Temp table
temp_table <- filter(frag_kegg_final, DE_origin == "temp", PValue < .05) |>  
  select(accession, 
         Description, 
         proteinId, 
         K_no, 
         logFC, 
         logCPM, 
         PValue)

temp_table[is.na(temp_table)] <- "-"

write.csv(temp_table, file =  here("Fig_S5-S7/CMA_NewPhytologist_temp_table_TMT_17052025.csv"))

#Int table
int_table <- filter(frag_kegg_final, DE_origin == "int", PValue < .05) |> 
  select(accession, 
         Description, 
         proteinId, 
         K_no, 
         logFC, 
         logCPM, 
         PValue)

int_table[is.na(int_table)] <- "-"

write.csv(int_table,  here("Fig_S5-S7/CMA_NewPhytologist_int_table_TMT_17052025.csv"))



# What proportion of proteins have no annotation? About half - 0.5330882
sum(!complete.cases(frag_kegg_final$K_no))/nrow(frag_kegg_final)

# Load in annotation text at KEGG levels
kegg_annotation <- read.csv(here("Fig_3/Fig_3_raw_data/KEGG2.csv")) |> 
  # Rename column to match
  rename("K_no" = "KO")


# Match by KO's
kegg_annotated <- left_join(frag_kegg_final, kegg_annotation, by = "K_no") |> 
  select(accession, 
         proteinId, 
         Description, 
         logFC, 
         logCPM, 
         PValue, 
         DE_origin, 
         de_dir, 
         K_no, 
         A, 
         B, 
         C, 
         D, 
         E, 
         FF,
         G, 
         H) |> 
  # Change NA's to unknown
  mutate(A = ifelse(is.na(A), 
                    'Unknown', A))

write.csv(kegg_annotated, here("Fig_3/Fig_3_output_tables/kegg_annotated_17052025.csv"))


# Barchart for DE'd proteins showing dist of functional annotations -------------------------------------------------

# Calculate counts
counts_df <- kegg_annotated |>  
  dplyr::count(DE_origin, 
               A, 
               B) |> 
  mutate(B = ifelse(is.na(B), 
                'No annotation available',
                B)) 

counts_df_B12 <- counts_df |> 
  filter(DE_origin == "B12")


# Get list of unique KEGG level B annotations from the counts
b_annotations <- unique(counts_df$B)


# Add missing rows 
missing_B12 <- setdiff(b_annotations,
                       counts_df_B12$B)

to_add_B12 <- filter(counts_df, B %in% missing_B12) |> 
  distinct(B, .keep_all = TRUE) |> 
  mutate(DE_origin = "B12",
         n = 0)
  


counts_df_B12 <- rbind(counts_df_B12, 
                       to_add_B12)


# Add a row with total of all counts
sum_row <- c("B12", NA ,"Total", sum(counts_df_B12$n))

counts_df_B12 <- counts_df_B12 |> 
  rbind(sum_row) |> 
  mutate(n = as.numeric(n)) |> 
  arrange(A, B) 


# Load color pallette
cols <- wes_palette("Royal2")


# Create list of lables to exlude 0's from geom_text argument
num_labels_b12 <- counts_df_B12$n
num_labels_b12[num_labels_b12 == "0"] <- " "
counts_df_B12$n[22:23] <- 0

# lock in factor level order
counts_df_B12$B <- factor(counts_df_B12$B, levels = counts_df_B12$B)


# Plot B12 Barchart ~ -----------------------------------------------------

# FIXME
# Wishlist - multipanel plot with labels
# zebra highlighting in background for functional groups 
# Fix labels in mulitpanel for subscripts on B12 etc


B12_plot <- ggplot(counts_df_B12, 
                   aes(y = B, 
                       x = n, 
                       fill = A)) + 
  geom_bar(position="dodge", 
           stat="identity") +
  theme_classic() +
  ylab("Functional Annotation") +
  xlab(NULL) +
  scale_y_discrete(limits = rev(levels(counts_df_B12$B))) +
  scale_fill_manual("Pathway", values = cols) +
  geom_text(aes(label= num_labels_b12),
            position=position_dodge(width=0.9),
            hjust=-.45, 
            size = 8) +
  theme(axis.text.y = element_text(size = 33),
        axis.title.y = element_text(size = 38, vjust = 5)) + 
  theme(text = element_text(size = 20), legend.position = "none", 
        #axis.text.y =element_blank()
        ) +
  theme(plot.margin = unit(c(3, 0, 3, 3), "cm")) +
  xlim(0,10) 



# Plot for Temp barchart ~ -------------------------------------------------

# Filter out for temp obs
counts_df_temp <- filter(counts_df, DE_origin == "temp")

# Add missing rows 
missing_temp <- setdiff(b_annotations,counts_df_temp$B)
to_add_temp <- filter(counts_df, B %in% missing_temp) |>  
  distinct(B, .keep_all = TRUE) |> 
  mutate(DE_origin = "temp", 
         n = 0)

counts_df_temp <- rbind(counts_df_temp, to_add_temp)


# Add a row with total of all counts
sum_row <- c("temp", NA ,"Total", sum(counts_df_temp$n))
counts_df_temp <- rbind(counts_df_temp, sum_row) |> 
  mutate(n = as.numeric(n))



# # Reorder obs 
counts_df_temp <- counts_df_temp |>  
  arrange(A, B) |> 
  mutate(B = factor(B, levels = B))



# Load color pallette
cols <- wes_palette("Royal2")


# Create list of labels to exlude 0's from geom_text argument
num_labels_temp <- counts_df_temp$n
num_labels_temp[num_labels_temp == "0"] <- " "


counts_df_temp$n[22:23] <- 0


# Plot temp !!! :-)
temp_plot <- ggplot(counts_df_temp, aes(y = B, 
                                        x = n, 
                                        fill = A)) +
  geom_bar(position = "dodge", 
           stat = "identity") +
  theme_classic() +
  ylab(NULL) +
  xlab(NULL) +
  scale_y_discrete(limits = rev(levels(counts_df_temp$B))) +
  scale_fill_manual("Pathway", values = cols) +
  geom_text(
    aes(label = num_labels_temp),
    position = position_dodge(width = 0.9),
    hjust = -.45,
    size = 8) +
  theme(plot.title = element_text(size = 20)) +
  theme(text = element_text(size = 20), axis.text.y = element_blank()) +
  xlim(0, 10) +
  theme(plot.margin = unit(c(3, 0, 3, 0), "cm"))

# FIXME
# Save size natively 
# ggsave("temp_annotation.pdf", width = 24, height = 18, units = "in")

# INT Plot -------------------------------------------------


# Filter out for int obs
counts_df_int <- filter(counts_df, DE_origin == "int")

# Add missing rows 
missing_int <- setdiff(b_annotations,counts_df_int$B)
to_add_int <- filter(counts_df, B %in% missing_int)
to_add_int <- to_add_int |> distinct(B, .keep_all = TRUE) |> 
  mutate(DE_origin = "int",
         n = 0)

counts_df_int <- rbind(counts_df_int, to_add_int)


# Add a row with total of all counts
sum_row <- c("int", NA ,"Total", sum(counts_df_int$n))
counts_df_int <- rbind(counts_df_int, sum_row) |> 
  mutate(n = as.numeric(n))



# # Reorder obs 
counts_df_int <- counts_df_int |> 
  arrange(A, B) |> 
  # Wrap text descriptions
  mutate(B = factor(B,
                    levels = B))





# Load color palette
cols <- wes_palette("Royal2")


# Create list of labels to exlude 0's from geom_text argument
num_labels_int <- counts_df_int$n
num_labels_int[num_labels_int == "0"] <- " "


counts_df_int$n[22:23] <- 0





# Plot Int barchart ~ -----------------------------------------------------

int_plot <- ggplot(counts_df_int, aes(y = B, x = n, fill = A)) + 
  geom_bar(position="dodge", stat="identity") +
  theme_classic() +
  ylab(NULL) +
  xlab(NULL) +
  xlim(0,30) +
  scale_y_discrete(limits = rev(levels(counts_df_int$B))) +
  scale_fill_manual("Pathway", values = cols) +
  geom_text(aes(label= num_labels_int), position=position_dodge(width=0.9), 
            hjust=-.45, size = 8) +
  theme(plot.title = element_text(size = 20)) + 
  theme(text = element_text(size=20), 
        axis.text.y =element_blank()
        ) +
  theme(plot.margin = unit(c(3, 0, 3, 0), "cm"), legend.position = "none") 



  
#Make a plot with just x-axis text
annotation_barplot <- ggarrange(B12_plot,
                                temp_plot,
                                int_plot, 
                                ncol = 3, 
                                nrow = 1,
                                common.legend = TRUE) + 
  rremove("ylab") + 
  rremove("xlab")


# Multipanel plot ---------------------------------------------------------

# Annotations for plot labels 
ann1 <- geom_text(aes(x = 0, 
                      y = 0, 
                      label = "-B[12]"),
                  parse = TRUE, 
                  size = 6)

ann2 <- geom_text(aes(x = 0, 
                      y = 0, 
                  label = "+12 °C)"),
                  parse = TRUE, 
                  size = 6)

ann3 <- geom_text(aes(x =0,
                      y = 0, 
                      label = "-B[12] and +12 °C)"),
                  parse = TRUE, 
                  size = 6)

# Arrange plots
# FIXME labels not rendering 
# Caused by error in `.f()`:
ggarrange(
  B12_plot,
  temp_plot,
  int_plot,
  ann1, # cannot convert to grob?? >:(
  ann2, 
  ann3,
  nrow = 2, 
  ncol = 3,
    # label.x = c(.75, .5, .5),
    # label.y = .975,
    common.legend = TRUE,
    widths = c(3.7, 1, 3),
    legend = "bottom"
  ) 



ggsave(here("Fig_3/Fig_3_plots/annotation_barchart_17052025.pdf"), 
       width = 35, 
       height = 18, 
       units = "in", 
       limitsize = FALSE)
