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

