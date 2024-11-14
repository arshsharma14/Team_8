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
abundance_file <- "../qiime/Soil/picrust2/pathway_abundance.tsv"
abundance_data <- fread(abundance_file, sep = "\t", header = TRUE, strip.white = TRUE)
abundance_data  =as.data.frame(abundance_data)

setnames(abundance_data, old = "#OTU ID", new = "pathway")

#Import your metadata file, no need to filter yet
metadata <- read_delim("../wrangling/soil_metadata.tsv")

setnames(metadata, old = "#SampleID", new = "sample-id")
#Example Looking at subject number
#If you have multiple variants, filter your metadata to include only 2 at a time

# Remove rows where cn_category is "Very Low", "Very High", or NULL (NA) in metadata
metadata <- metadata[!(metadata$cn_category %in% c("Very Low", "Very High")) & !is.na(metadata$cn_category), ]

#Remove NAs for your column of interest in this case subject
metadata = metadata[!is.na(metadata$cn_category),]

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
abundance_desc = abundance_desc[,-c(574:ncol(abundance_desc))] 

# Generate pathway PCA plot
PCA <- pathway_pca(abundance = abundance_data_filtered %>% column_to_rownames("pathway"), metadata = metadata, group = "cn_category")

ggsave("soil_pca.png"
       , PCA
       , height=8, width =12)
