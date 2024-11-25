# Script for functional annotation of DE proteins from TMT experiment 


# Setup :-) ---------------------------------------------------------------

# Set working dir 
setwd("~/Bertrand Lab Dropbox/Bertrand Lab shared workspace/Catalina/Summer_2022/0_Fracy_TMT_Manuscript/Fracy_TMT_FigScripts/Fig_9")


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

# Installing KEGGREST
# BiocManager::install("KEGGREST")

# 4,B12 vs 4,noB12  (B12 response)
B12_hits <- read.csv("hits_4_noB12_20102021.csv") |>
  mutate(Description = str_remove(Description, fixed(pattern = "[Fragilariopsis cylindrus CCMP1102]")),
         DE_origin = "B12")



# 4,B12 vs 12,B12 (temp response)
temp_hits <- read.csv("hits_12_B12_20102021.csv") |>
  mutate(Description = str_remove(Description, fixed(pattern = "[Fragilariopsis cylindrus CCMP1102]")),
         DE_origin = "temp")

# 4,B12 vs 12,noB12 (interaction repsponse between B12 and temp)
int_hits <- read.csv("hits_12_noB12_20102021.csv") |>
mutate(Description = str_remove(Description, fixed(pattern = "[Fragilariopsis cylindrus CCMP1102]")),
       DE_origin = "int")


# Merge hit lists into 1, only including significantly DE'd proteins (Pvalues < .05)
hits_list <- rbind(B12_hits, temp_hits, int_hits) |> filter(PValue < .05) |>
  # Make a column that indicates DE direction (up or downregulated; positive or negative fold change)
  mutate(de_dir = fifelse(logFC < 0,
                          -1,
                          1))

# Add KEGG Annotation Information -----------------------------------------

# Match in protein ID's from NCBI (https://www.ncbi.nlm.nih.gov/genome/browse/#!/proteins/11246/283952%7CFragilariopsis%20cylindrus%20CCMP1102/)
prot_info <- read.csv("prot_table.csv") |>
  rename("accession" = "Protein.product")

# Protein ID's for chloroplast from (https://www.ncbi.nlm.nih.gov/genome/browse/#!/proteins/11246/748136%7CFragilariopsis%20cylindrus/chloroplast/)
# chl_prot_info <- read.csv("chl_prot_table.csv")

# Match accessions to protein ID's
prot_data_ids <- left_join(hits_list, prot_info, by = "accession") |>
  
  # Make new column with just protein ID's (6 digit code at end of locus.tag)
  mutate(proteinId = gsub(".*_", "", Locus.tag))


# Match ID's to functional information from KEGG
# Use KEGGREST package to access KEGG REST API 
# https://www.bioconductor.org/packages/release/bioc/vignettes/KEGGREST/inst/doc/KEGGREST-vignette.html
KEGGlist_frag <- as.data.frame(keggList("fcy"))
frag_kegg_nos_raw <- read.csv("frag_kegg_no.csv")

frag_kegg_nos <- colsplit(frag_kegg_nos_raw$Name," ",c("id","kegg_id"))


thing <- colsplit(frag_kegg_nos$kegg_id," ", c("thing","description"))
database <- cbind(frag_kegg_nos, thing)

# Delete row if no KO number is available 
frag_kegg <- database[!grepl("no", database$thing),]

# Remove second unneeded column 
frag_kegg <- frag_kegg[,-c(2)]

#Rename column with K numbers
colnames(frag_kegg)[2] <- "K_no" 

# Isolate protein Id's from column
prot_id_df <- colsplit(frag_kegg$id,"_",c("fracy","prot_id"))

# Add back into df
frag_kegg_database <- cbind(prot_id_df$prot_id, frag_kegg)

# Remove unneeded 2nd column and rename to produce a dataframe matching protein id's and kegg numbers
colnames(frag_kegg_database)[1] <- "proteinId"
frag_kegg_database <- frag_kegg_database[,-c(2,4)]

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

write.csv(B12_table, 
          file = "B12_table_TMT_29102021.csv")


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
write.csv(temp_table, file = "temp_table_TMT_29102021.csv")

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

write.csv(int_table, file = "int_table_TMT_29102021.csv")



# What proportion of proteins have no annotation? About half
sum(!complete.cases(frag_kegg_final$K_no))/nrow(frag_kegg_final)

# Load in annotation text at KEGG levels
kegg_annotation <- read.csv("KEGG2.csv")

# Change column name
colnames(kegg_annotation)[which(names(kegg_annotation) == "KO")] <- "K_no"

# Match by KO's
frag_kegg_annotated <- left_join(frag_kegg_final, kegg_annotation, by = "K_no")

# Select only useful cols
kegg_annotated <- frag_kegg_annotated |> 
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
         H)


# Plot for B12 ------------------------------------------------------------

# Create a table with annotated proteins
B12_de_annotated <- kegg_annotated |>  filter(DE_origin == "B12") |> 
  mutate(A = ifelse(is.na(A), 
                    'Unknown', A))



THING <- B12_de_annotated %>% group_by(accession, A) %>% dplyr::summarise(logFC= first(logFC), logCPM = first(logCPM))


b12_accessions <- unique(B12_de_annotated$accession)

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
b12_labels_df <- as.data.frame(cbind(b12_accessions, b12_names))

colnames(b12_labels_df) <- c("accession", "label")

b12_de_annotated_labeled <- THING |> 
  left_join(b12_labels_df, by = "accession")

# add chld 
g1 = subset(int_de_annotated_labeled, accession == "OEU10229.1")

g2 = subset(int_de_annotated_labeled, accession == "OEU11144.1")


b12_scatter <- ggplot(b12_de_annotated_labeled, aes(x = logCPM, y = logFC, color =
                                               A)) +
  geom_point(size = 5, alpha = .50) +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  xlab(expression('Log'[2] * ' Abundance')) +
  ylab(expression('Log'[2] * ' Fold Change')) +
  scale_color_manual(values = wes_palette("Royal2"), name = "Pathway") +
  geom_text_repel(data = g1, 
                  label = "ChlD", 
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

temp_de_annotated <- kegg_annotated %>% filter(DE_origin == "temp")
temp_de_annotated$A <- ifelse(is.na(temp_de_annotated$A), 
                              'Unknown', temp_de_annotated$A)

temp_label <- c("Unknown Protein (ID: 271832)", "Unknown Protein (ID: 235337)")

temp_scatter <- ggplot(temp_de_annotated, aes(x = logCPM, y = logFC, color =
                                                A)) +
  geom_point(size = 5, alpha = .50) +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  xlab(expression('Log'[2] * ' Abundance')) +
  ylab(expression('Log'[2] * ' Fold Change')) +
  scale_color_manual(values = wes_palette("Royal2"), name = "Pathway") +
  geom_text_repel(
    aes(label = ifelse(logFC > 2.5 & logCPM > 8, temp_label, '')),
    
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
int_de_annotated <- kegg_annotated %>% filter(DE_origin == "int")
int_de_annotated$A <- ifelse(is.na(int_de_annotated$A), 
                             'Unknown', int_de_annotated$A)

# Create table with shorter alt names for important proteins (logFC >/< 2.5; logCPM > 12)
accession <- c("OEU18748.1", "OEU11144.1", "OEU08040.1", "OEU11214.1", "OEU19288.1", "OEU15214.1", "OEU08856.1", "OEU22459.1", "OEU13453.1", "OEU15467.1", "OEU16276.1", "OEU11005.1")

names <- c("Unknown Protein (ID: 182880)", "MetE", "Unknown Protein (ID: 271832)", "CBA1", "Unknown Protein (ID: 235337)", "Nucleoside Triphosphate Hydrolase", "EngA", "NADH oxireductase", "Dihydrolipoamide acetyltransferase", "Triosephosphate isomerase", "Asparagine Synthase", "Band7 Domain Containing Protein" )

labels_df <- as.data.frame(cbind(accession, names))


int_de_annotated_labeled <- int_de_annotated %>% left_join(labels_df, by = "accession")

# Remove second MetE enrty to label
int_de_annotated_labeled <- int_de_annotated_labeled[-c(3),]

# Add ChlD for discussion purposes
g1 <- subset(int_de_annotated_labeled, accession == "OEU10229.1")

# Plotty plotty 
int_scatter <- ggplot(int_de_annotated_labeled, aes(x=logCPM, y=logFC, color=A)) +
  geom_point(size = 5, alpha = .50) +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  xlab(expression('Log'[2]*' Abundance')) +
  ylab(expression('Log'[2]*' Fold Change')) +
  scale_color_manual(values = wes_palette("Royal2"), name = "Pathway")+
  geom_text_repel(aes(label=ifelse(logFC < -2.5 & logCPM > 8,as.character(names),'')),
                  #hjust=0,
                  #   vjust=0, 
                  show.legend =  FALSE, 
                  nudge_y = -1, 
                  nudge_x = 1,
                  size = 5) +
  geom_text_repel(aes(label=ifelse(logFC > 2.5& logCPM > 8,
                                   as.character(names),'')),
                  #        hjust=0,
                  #vjust=0, 
                  show.legend =  FALSE, 
                  nudge_y = 1, 
                  nudge_x = 1.75,
                  size = 5) +  
  geom_text_repel(data = g1, 
                  label = "ChlD", 
                  vjust = -1.5, 
                  hjust = -.75,
                  show.legend = FALSE, 
                  size = 5, 
                  nudge_x = 1.5) +
  xlim(7,20) +
  ylim(-7, 7)+
  ggtitle(expression("-B"[12]*" and +12° C")) +
  theme(plot.title = element_text(hjust=0.5), 
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )


plot <- ggarrange(b12_scatter, temp_scatter, int_scatter, ncol = 3, nrow = 1, common.legend = TRUE, legend = "bottom")

ggsave(plot, filename = "scatter_all_2_24112024.png",  bg = "transparent", width = 16, height= 7, units = "in")


# barchart for DE'd proteins showing dist of functional annotations -------------------------------------------------


# Calculate counts
counts_df <- frag_kegg_annotated %>% dplyr::count(DE_origin, A, B)

# Get list of unique KEGG level B annotations from the counts
b_annotations <- unique(counts_df$B)


# Change NA's to text
counts_df$B <- ifelse(is.na(counts_df$B), 
                      'No annotation available', counts_df$B)
counts_df$A <- ifelse(is.na(counts_df$A), 
                      'Unknown', counts_df$A)


# Filter out for B12 obs
counts_df_B12 <- filter(counts_df, DE_origin == "B12")

# Add missing rows 
missing_B12 <- setdiff(b_annotations,counts_df_B12$B)
to_add_B12 <- filter(counts_df, B %in% missing_B12)
to_add_B12<- to_add_B12 %>% distinct(B, .keep_all = TRUE)
to_add_B12$DE_origin <- "B12"
to_add_B12$n <- 0
counts_df_B12 <- rbind(counts_df_B12, to_add_B12)


# Add a row with total of all counts
sum_row <- c("B12", NA ,"Total", sum(counts_df_B12$n))
counts_df_B12 <- rbind(counts_df_B12, sum_row)
counts_df_B12$n <- as.numeric(counts_df_B12$n )


# Reorder obs 
counts_df_B12 <- counts_df_B12[
   with(counts_df_B12, order(A, B)),]

# Make text wrapping
counts_df_B12$B <- stringr::str_wrap(counts_df_B12$B, 35)

# lock in factor level order
counts_df_B12$B <- factor(counts_df_B12$B, levels = counts_df_B12$B)

# Load color pallette
cols <- wes_palette("Royal2")


# Create list of lables to exlude 0's from geom_text argument
num_labels_b12 <- counts_df_B12$n
num_labels_b12[num_labels_b12 == "0"] <- " "


counts_df_B12$n[22:23] <- 0

# Plot B12 !!! :-)
B12_plot <- ggplot(counts_df_B12, aes(y = B, x = n, fill = A)) + 
  geom_bar(position="dodge", stat="identity") +
  theme_classic() +
  ylab("Functional Annotation") +
  xlab(NULL) +
  scale_y_discrete(limits = rev(levels(counts_df_B12$B))) +
  scale_fill_manual("Pathway", values = cols) +
  xlim(0,30) +
  geom_text(aes(label= num_labels_b12), position=position_dodge(width=0.9), 
            hjust=-.45, size = 8) +
  theme(plot.title = element_text(size = 20)) + 
  theme(text = element_text(size=40), legend.position = "none", 
        #axis.text.y =element_blank()
        ) +
  theme(plot.margin = unit(c(3, 0, 3, 3), "cm"))

ggsave("b12_annotation.pdf", width = 24, height = 18, units = "in")

# Temp Plot -------------------------------------------------

# Filter out for temp obs
counts_df_temp <- filter(counts_df, DE_origin == "temp")

# Add missing rows 
missing_temp <- setdiff(b_annotations,counts_df_temp$B)
to_add_temp <- filter(counts_df, B %in% missing_temp)
to_add_temp<- to_add_temp %>% distinct(B, .keep_all = TRUE)
to_add_temp$DE_origin <- "temp"
to_add_temp$n <- 0
counts_df_temp <- rbind(counts_df_temp, to_add_temp)


# Add a row with total of all counts
sum_row <- c("temp", NA ,"Total", sum(counts_df_temp$n))
counts_df_temp <- rbind(counts_df_temp, sum_row)
counts_df_temp$n <- as.numeric(counts_df_temp$n )


# # Reorder obs 
counts_df_temp <- counts_df_temp %>% arrange(A, B)

# Make text wrapping
counts_df_temp$B <- stringr::str_wrap(counts_df_temp$B, 35)



# Load color pallette
cols <- wes_palette("Royal2")


# Create list of labels to exlude 0's from geom_text argument
num_labels_temp <- counts_df_temp$n
num_labels_temp[num_labels_temp == "0"] <- " "


counts_df_temp$n[22:23] <- 0

# lock in factor level order
counts_df_temp$B <- factor(counts_df_temp$B, levels = counts_df_temp$B)

# Plot temp !!! :-)
temp_plot <- ggplot(counts_df_temp, aes(y = B, x = n, fill = A)) + 
  geom_bar(position="dodge", stat="identity") +
  theme_classic() +
  ylab(NULL) +
  xlab(NULL) +
  scale_y_discrete(limits = rev(levels(counts_df_temp$B))) +
  scale_fill_manual("Pathway", values = cols) +
  geom_text(aes(label= num_labels_temp), position=position_dodge(width=0.9), 
            hjust=-.45, size = 8) +
  theme(plot.title = element_text(size = 20)) + 
  theme(text = element_text(size=40),
       axis.text.y =element_blank()
  ) +
  xlim(0,30) +
  theme(plot.margin = unit(c(3, 0, 3, 0), "cm"))


ggsave("temp_annotation.pdf", width = 24, height = 18, units = "in")

# INT Plot -------------------------------------------------


# Filter out for int obs
counts_df_int <- filter(counts_df, DE_origin == "int")

# Add missing rows 
missing_int <- setdiff(b_annotations,counts_df_int$B)
to_add_int <- filter(counts_df, B %in% missing_int)
to_add_int <- to_add_int %>% distinct(B, .keep_all = TRUE)
to_add_int$DE_origin <- "int"
to_add_int$n <- 0
counts_df_int <- rbind(counts_df_int, to_add_int)


# Add a row with total of all counts
sum_row <- c("int", NA ,"Total", sum(counts_df_int$n))
counts_df_int <- rbind(counts_df_int, sum_row)
counts_df_int$n <- as.numeric(counts_df_int$n )


# # Reorder obs 
counts_df_int <- counts_df_int %>% arrange(A, B)

# Make text wrapping
counts_df_int$B <- stringr::str_wrap(counts_df_int$B, 35)



# Load color palette
cols <- wes_palette("Royal2")


# Create list of labels to exlude 0's from geom_text argument
num_labels_int <- counts_df_int$n
num_labels_int[num_labels_int == "0"] <- " "


counts_df_int$n[22:23] <- 0

# # lock in factor level order
 counts_df_int$B <- factor(counts_df_int$B, levels = counts_df_int$B)


# Plot int !!! :-)
int_plot <- ggplot(counts_df_int, aes(y = B, x = n, fill = A)) + 
  geom_bar(position="dodge", stat="identity") +
  theme_classic() +
  ylab(NULL) +
  xlab(NULL) +
  scale_y_discrete(limits = rev(levels(counts_df_int$B))) +
  scale_fill_manual("Pathway", values = cols) +
  geom_text(aes(label= num_labels_int), position=position_dodge(width=0.9), 
            hjust=-.45, size = 8) +
  theme(plot.title = element_text(size = 20)) + 
  theme(text = element_text(size=40), 
     #   axis.text.y =element_blank()
        ) +
  theme(plot.margin = unit(c(3, 0, 3, 0), "cm"), legend.position = "none") 


ggsave("int_annotation.pdf", width = 24, height = 18, units = "in")
  
#Make a plot with just x-axis text
# Dummy values to plot 
dummy <- data.frame(levels(counts_df_int$B))





ggarrange(B12_plot, temp_plot, int_plot, ncol = 3, nrow = 1, common.legend = TRUE)
