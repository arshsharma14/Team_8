# open libraries
library(tidyverse)
library(readxl)

# open soil metadata
soil_meta <- read_xlsx("soil_metadata.xlsx")

# select for needed columns
soil_select <- select(soil_meta, `#SampleID`, `Sample Accession`, `Sample Alias`, Ecozone, Region, Site, 
                      `Environmental Source`, `Phylogenetic Marker`, `Target Gene Subfragment`, `Collection Date`, 
                      `Sampling Depth`, Elevation, `Mean Annual Temperature (Celsius)`, `Total Carbon`, `Total Nitrogen`, pH, `CN Ratio`)

soil_select <- rename(soil_select, `elevation_m` = `Elevation`, `sampling_depth_m` = `Sampling Depth`, 
                      `mean_annual_temp_C` = `Mean Annual Temperature (Celsius)`, `ph` = `pH`, `cn_ratio` = 
                        `CN Ratio`, `total_carbon` = `Total Carbon`, `total_nitrogen` = `Total Nitrogen`)                 

soil_select$ph <- as.numeric(soil_select$ph)


soil_bins <- c(-Inf, 3, 5, 6.5, 7.5, 9, 11, Inf)

soil_bin_labels <- c("Strongly Acidic", "Moderately Acidic", "Slightly Acidic", 
                    "Neutral", "Slightly Basic", "Moderately Basic", "Strongly Basic")

soil_select$ph[soil_select$ph == 0.0] <- NA # Some of the 0.0 values in ph column are supposed to be NA
sum(is.na(soil_select$ph))
sum(!is.na(soil_select$ph))

soil_select$ph_category <- cut(soil_select$ph, # make new column "ph_category"
                               breaks = soil_bins,        # cut function: divides numeric into intervals and labels them according to assigned labels
                               labels = soil_bin_labels, 
                               right = FALSE)