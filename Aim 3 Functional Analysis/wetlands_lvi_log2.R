#### Load packages ####
# Load all necessary libraries
library(readr)
library(data.table)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(DESeq2)
library(ggh4x)


#### Import files and preparing tables ####
#Importing the pathway PICrsut2
abundance_file <- "../qiime/Wetlands/picrust2/pathway_abundance.tsv"
abundance_data <- fread(abundance_file, sep = "\t", header = TRUE, strip.white = TRUE)
abundance_data  =as.data.frame(abundance_data)

setnames(abundance_data, old = "#OTU ID", new = "pathway")

#Import your metadata file, no need to filter yet
metadata <- read_delim("../wrangling/wetlands_metadata.tsv")

setnames(metadata, old = "#SampleID", new = "sample-id")
#Example Looking at subject number
#If you have multiple variants, filter your metadata to include only 2 at a time

#Remove NAs for your column of interest in this case subject
metadata = metadata[!is.na(metadata$cn_category),]

# Filter for only "Low" and "Intermediate" in cn_category
metadata <- metadata[metadata$cn_category %in% c("Low", "Intermediate"), ]

#Filtering the abundance table to only include samples that are in the filtered metadata
sample_names = metadata$'sample-id'
sample_names = append(sample_names, "pathway")
abundance_data_filtered = abundance_data[, colnames(abundance_data) %in% sample_names] #This step is the actual filtering

#Removing individuals with no data that caused a problem for pathways_daa()
abundance_data_filtered =  abundance_data_filtered[, colSums(abundance_data_filtered != 0) > 0]

#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(abundance_data_filtered) = NULL

#verify samples in metadata match samples in abundance_data
abun_samples = rownames(t(abundance_data_filtered[,-1])) #Getting a list of the sample names in the newly filtered abundance data
metadata = metadata[metadata$`sample-id` %in% abun_samples,] #making sure the filtered metadata only includes these samples

#### DESEq ####
#Perform pathway DAA using DESEQ2 method
abundance_daa_results_df <- pathway_daa(abundance = abundance_data_filtered %>% column_to_rownames("pathway"), 
                                        metadata = metadata, group = "cn_category", daa_method = "DESeq2")

# Annotate MetaCyc pathway so they are more descriptive
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", 
                                                       daa_results_df = abundance_daa_results_df, ko_to_kegg = FALSE)

# Filter p-values to only significant ones
feature_with_p_0.05 <- abundance_daa_results_df %>% filter(p_values < 0.05)

#Changing the pathway column to description for the results 
feature_desc = inner_join(feature_with_p_0.05,metacyc_daa_annotated_results_df, by = "feature")
feature_desc$feature = feature_desc$description
feature_desc = feature_desc[,c(1:7)]
colnames(feature_desc) = colnames(feature_with_p_0.05)

#Changing the pathway column to description for the abundance table
abundance = abundance_data_filtered %>% filter(pathway %in% feature_with_p_0.05$feature)
colnames(abundance)[1] = "feature"
abundance_desc = inner_join(abundance,metacyc_daa_annotated_results_df, by = "feature")
abundance_desc$feature = abundance_desc$description
#this line will change for each dataset. 34 represents the number of samples in the filtered abundance table
abundance_desc = abundance_desc[,-c(94:ncol(abundance_desc))] 

# Filter tables for category vs category - abundance_desc ; abundance_data_filtered ; metacyc_daa_annotated_results_df

# Generate a heatmap
heatmap <- pathway_heatmap(abundance = abundance_desc %>% column_to_rownames("feature"), metadata = metadata, group = "cn_category")

# Generating a bar plot representing log2FC from the custom deseq2 function

# Go to the Deseq2 function script and update the metadata category of interest

# Lead the function in
source("DESeq2_function.R")

# Run the function on your own data
res =  DEseq2_function(abundance_data_filtered, metadata, "cn_category")
res$feature =rownames(res)
res_desc = inner_join(res,metacyc_daa_annotated_results_df, by = "feature")
res_desc = res_desc[, -c(8:13)]

# Filter to only include significant pathways
sig_res = res_desc %>%
  filter(pvalue < 0.05)

# Filter log2FoldChange values to keep only those > 1 or < -1
sig_res_1 <- sig_res[sig_res$log2FoldChange > 1 | sig_res$log2FoldChange < -1, ]

sig_res_1 <- sig_res_1[order(sig_res_1$log2FoldChange),]
log2_fold_change_1 <- ggplot(data = sig_res_1, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(
    title = "Wetlands Low vs Intermediate", 
    x = "Log2 Fold Change", 
    y = "Pathways",
    fill = expression(bold("P Value"))  
  ) +
  theme(
    plot.title = element_text(size = 22, face = "bold"),  
    axis.text = element_text(size = 14),                  
    axis.title = element_text(size = 18, face = "bold")
  )
log2_fold_change_1

ggsave("wetlands_lvi_log2_1.png"
       , log2_fold_change_1
       , height=14, width =12)

# Count the number of positive log2FoldChange values
positive_count <- sum(sig_res_1$log2FoldChange > 0)

# Count the number of negative log2FoldChange values
negative_count <- sum(sig_res_1$log2FoldChange < 0)

# Create a new data frame 
count <- data.frame(
  log2FoldChange = c("Positive", "Negative"),
  Count = c(positive_count, negative_count)
)
count

