# open libraries
library(tidyverse)
library(readxl)

# open soil metadata
meta <- read_xlsx("soil_metadata.xlsx")

# select for needed columns
meta_select <- select(meta, `#SampleID`, `Sample Accession`, `Sample Alias`, Ecozone, Region, Site, `Environmental Source`, `Phylogenetic Marker`, `Target Gene Subfragment`, `Collection Date`, `Sampling Depth`, Elevation, `Mean Annual Temperature (Celsius)`, `Total Carbon`, `Total Nitrogen`, pH, `CN Ratio`)

meta_select <- rename(meta_select, `elevation_m` = `Elevation`, `sampling_depth_m` = `Sampling Depth`, `mean_annual_temp_C` = `Mean Annual Temperature (Celsius)`, `ph` = `pH`, `cn_ratio` = `CN Ratio`, `total_carbon` = `Total Carbon`, `total_nitrogen` = `Total Nitrogen`)                 
