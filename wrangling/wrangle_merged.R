# open libraries
library(tidyverse)

metadata <- read_tsv("merged-metadata.tsv", )
metadata <- metadata[-1,]

# Remove rows where cn_category is "Very Low", "Very High", or NULL (NA) in metadata
metadata <- metadata[!(metadata$cn_category %in% c("Very Low", "Very High")) & !is.na(metadata$cn_category), ]

category_counts_cn <- table(metadata$cn_category)
category_counts_df_cn <- as.data.frame(category_counts_cn)

# Replace NA values in Environmental Source with "Wetlands"
metadata$"Environmental Source"[is.na(metadata$"Environmental Source")] <- "Wetlands"

# Create the new column 'cn_category_per_dataset' by combining 'cn_category' and 'Environmental Source'
metadata$cn_category_per_dataset <- ifelse(
  metadata$"Environmental Source" == "Forest Soil",
  paste(metadata$cn_category, "Soil", sep = "-"),
  ifelse(metadata$"Environmental Source" == "Wetlands",
         paste(metadata$cn_category, "Wetlands", sep = "-"),
         NA)
)

category_counts_cn_pdf <- table(metadata$cn_category_per_dataset)
category_counts_df_cn_pdf <- as.data.frame(category_counts_cn_pdf)

write.table(metadata, file = "merged_wrangled_metadata.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
