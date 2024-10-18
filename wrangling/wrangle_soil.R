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

category_counts_soil <- table(soil_select$ph_category)
category_counts_df_soil <- as.data.frame(category_counts_soil)

#Convert C:N ratio to numeric
soil_select$cn_ratio <- as.numeric(soil_select$cn_ratio)
range(soil_select$cn_ratio, na.rm = TRUE)
ggplot(soil_select, aes(x=cn_ratio)) +
  geom_histogram()

#looking at C:N ratio after removing NA and those over 50
soil_cn <- soil_select %>% filter(!is.na(cn_ratio) & cn_ratio < 50)
range(soil_cn$cn_ratio, na.rm = TRUE)
ggplot(soil_cn, aes(x=cn_ratio)) +
  geom_histogram()

#Adding bins for C:N ratio
cn_bins <- c(-Inf, 10, 20, 30, 50, Inf)

cn_bin_labels <- c("Very Low", "Low", "Intermediate", "High", "Very High")

soil_select$cn_category <- cut(soil_select$cn_ratio, # make new column "ph_category"
                               breaks = cn_bins,        # cut function: divides numeric into intervals and labels them according to assigned labels
                               labels = cn_bin_labels, 
                               right = FALSE)

ggplot(soil_select, aes(x=cn_category)) +
  geom_bar()

write.table(soil_select, file = "soil_metadata.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

max(soil_select$cn_ratio, na.rm = TRUE)