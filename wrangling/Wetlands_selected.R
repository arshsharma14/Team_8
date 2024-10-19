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

# Converting carbon and nitrogen content columns to numeric
wetlands_selected$total_carbon_percent <- as.numeric(wetlands_selected$total_carbon_percent)
wetlands_selected$total_nitrogen_percent <- as.numeric(wetlands_selected$total_nitrogen_percent)

# Calculating C:N ratio 
wetlands_selected$cn_ratio <- wetlands_selected$total_carbon_percent / wetlands_selected$total_nitrogen_percent
range(wetlands_selected$cn_ratio, na.rm = TRUE)

# Looking at C:N ratio distribution for whole metadata and for those without NA and > 50
ggplot(wetlands_selected, aes(x=cn_ratio)) +
  geom_histogram()

wetlands_cn <- wetlands_selected %>% filter(!is.na(cn_ratio) & cn_ratio < 50)
range(wetlands_cn$cn_ratio, na.rm = TRUE)
ggplot(wetlands_cn, aes(x=cn_ratio)) +
  geom_histogram()

# Adding bins for C:N ratio
cn_bins <- c(-Inf, 10, 20, 30, 50, Inf)

cn_bin_labels <- c("Very Low", "Low", "Intermediate", "High", "Very High")

wetlands_selected$cn_category <- cut(wetlands_selected$cn_ratio, # make new column "ph_category"
                               breaks = cn_bins,        # cut function: divides numeric into intervals and labels them according to assigned labels
                               labels = cn_bin_labels, 
                               right = FALSE)

ggplot(wetlands_selected, aes(x=cn_category)) +
  geom_bar()

category_counts_w_cn <- table(wetlands_selected$cn_category)
category_counts_df_w_cn <- as.data.frame(category_counts_w_cn)

write.table(wetlands_selected, file = "wetlands_metadata.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

max(wetlands_selected$cn_ratio, na.rm = TRUE)