install.packages("readxl")
library(readxl)
library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)

#load in wetlands metadata
meta <- "wetlands_metadata.xlsx"
wetlands <- read_excel(meta)

#selecting variables considered for project
wetlands_selected <- wetlands[, c("#SampleID", "env_feature", "env_material", "elevation_m", 
                               "total_organic_carb_percent","total_inorganic_carbon_percent",  
                               "total_carbon_percent", "soil_total_inorganic_nitrogen", 
                               "soil_nitrogen_units", "total_nitrogen_percent", "ph")]
str(wetlands_selected$ph)
wetlands_selected$ph <- as.numeric(wetlands_selected$ph)


min(wetlands_selected$ph, na.rm = TRUE)
max(wetlands_selected$ph, na.rm = TRUE)

bins <- c(-Inf, 3, 5, 6.5, 7.5, 9, 11, Inf) # Defines breaks/boundaries for each category

bin_labels <- c("Strongly Acidic", "Moderately Acidic", "Slightly Acidic", 
                "Neutral", "Slightly Basic", "Moderately Basic", "Strongly Basic")

wetlands_selected$ph_category <- cut(wetlands_selected$ph, # make new column "ph_category"
                                     breaks = bins,        # cut function: divides numeric into intervals and labels them according to assigned labels
                                     labels = bin_labels, 
                                     right = FALSE)

category_counts <- table(wetlands_selected$ph_category)
category_counts_df <- as.data.frame(category_counts)

