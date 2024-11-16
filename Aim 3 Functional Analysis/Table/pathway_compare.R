# load in libraries
library(dplyr)

# load in dataframes
load("soil_int_v_high_table.RData")
load("soil_low_v_high_table.RData")
load("soil_low_v_int_table.RData")

load("wetlands_int_v_high_table.RData")
load("wetlands_low_v_high_table.RData")
load("wetlands_low_v_int_table.RData")

# Find the maximum number of rows among the three dataframes
max_rows <- max(nrow(pathways1), nrow(pathways2), nrow(pathways3))

# Add rows of NA to each dataframe to make their number of rows equal to max_rows
pathways1_padded <- pathways1 %>% 
  add_row(.before = nrow(pathways1) + 1, .rows = max_rows - nrow(pathways1))

pathways2_padded <- pathways2 %>% 
  add_row(.before = nrow(pathways2) + 1, .rows = max_rows - nrow(pathways2))

pathways3_padded <- pathways3 %>% 
  add_row(.before = nrow(pathways3) + 1, .rows = max_rows - nrow(pathways3))

# Bind the columns together
soil_pathways <- bind_cols(pathways1_padded, pathways2_padded, pathways3_padded)

# Find the maximum number of rows among the three dataframes
max_rows <- max(nrow(pathways4), nrow(pathways5), nrow(pathways6))

# Add rows of NA to each dataframe to make their number of rows equal to max_rows
pathways4_padded <- pathways4 %>% 
  add_row(.before = nrow(pathways4) + 1, .rows = max_rows - nrow(pathways4))

pathways5_padded <- pathways5 %>% 
  add_row(.before = nrow(pathways5) + 1, .rows = max_rows - nrow(pathways5))

pathways6_padded <- pathways6 %>% 
  add_row(.before = nrow(pathways6) + 1, .rows = max_rows - nrow(pathways6))

# Bind the columns together
wetlands_pathways <- bind_cols(pathways4_padded, pathways5_padded, pathways6_padded)

# Find common pathways between the columns High_v_Int and High_v_Low
common_pathways <- intersect(soil_pathways$High_v_Int[!is.na(soil_pathways$High_v_Int)], 
                             soil_pathways$High_v_Low[!is.na(soil_pathways$High_v_Low)])

# Create a new column High with common pathways and pad the rest with NA
max_rows <- nrow(soil_pathways)
high_list <- c(common_pathways, rep(NA, max_rows - length(common_pathways)))

# Add the new column to the soil_pathways DataFrame
soil_pathways$High <- high_list

# Find common pathways between the columns Int_v_Low and Int_v_High
common_pathways_int <- intersect(soil_pathways$Int_v_Low[!is.na(soil_pathways$Int_v_Low)], 
                                 soil_pathways$Int_v_High[!is.na(soil_pathways$Int_v_High)])

# Create a new column Int with common pathways and pad the rest with NA
max_rows <- nrow(soil_pathways)
int_list <- c(common_pathways_int, rep(NA, max_rows - length(common_pathways_int)))

# Add the new column to the soil_pathways DataFrame
soil_pathways$Int <- int_list

# Find common pathways between the columns Low_v_Int and Low_v_High
common_pathways_low <- intersect(soil_pathways$Low_v_Int[!is.na(soil_pathways$Low_v_Int)], 
                                 soil_pathways$Low_v_High[!is.na(soil_pathways$Low_v_High)])

# Create a new column Low with common pathways and pad the rest with NA
max_rows <- nrow(soil_pathways)
low_list <- c(common_pathways_low, rep(NA, max_rows - length(common_pathways_low)))

# Add the new column to the soil_pathways DataFrame
soil_pathways$Low <- low_list

# Step 1: Create High Column (common pathways between High_v_Int and High_v_Low)
common_pathways_high <- intersect(wetlands_pathways$High_v_Int[!is.na(wetlands_pathways$High_v_Int)], 
                                  wetlands_pathways$High_v_Low[!is.na(wetlands_pathways$High_v_Low)])
max_rows <- nrow(wetlands_pathways)
high_list <- c(common_pathways_high, rep(NA, max_rows - length(common_pathways_high)))

# Step 2: Create Int Column (common pathways between Int_v_Low and Int_v_High)
common_pathways_int <- intersect(wetlands_pathways$Int_v_Low[!is.na(wetlands_pathways$Int_v_Low)], 
                                 wetlands_pathways$Int_v_High[!is.na(wetlands_pathways$Int_v_High)])
int_list <- c(common_pathways_int, rep(NA, max_rows - length(common_pathways_int)))

# Step 3: Create Low Column (common pathways between Low_v_Int and Low_v_High)
common_pathways_low <- intersect(wetlands_pathways$Low_v_Int[!is.na(wetlands_pathways$Low_v_Int)], 
                                 wetlands_pathways$Low_v_High[!is.na(wetlands_pathways$Low_v_High)])
low_list <- c(common_pathways_low, rep(NA, max_rows - length(common_pathways_low)))

# Add High, Int, and Low columns to df
wetlands_pathways$High = high_list
wetlands_pathways$Int = int_list
wetlands_pathways$Low = low_list