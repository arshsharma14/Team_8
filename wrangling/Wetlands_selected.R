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

