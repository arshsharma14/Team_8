# Core Microbiome Libraries #
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

#### Load data ####
load("wetlands_final.RData")

# Convert to relative abundance
mpt_RA <- transform_sample_counts(mpt_final, fun=function(x) x/sum(x))

# Filter dataset by antibiotic use
mpt_low <- subset_samples(mpt_RA, `cn_category`=="Low")
mpt_int <- subset_samples(mpt_RA, `cn_category`=="Intermediate")
mpt_high <- subset_samples(mpt_RA, `cn_category`=="High")

# What ASVs are found in more than 70% of samples in each C:N category?
# trying changing the prevalence to see what happens
low_ASVs <- core_members(mpt_low, detection=0, prevalence = 0.7)
int_ASVs <- core_members(mpt_int, detection=0, prevalence = 0.7)
high_ASVs <- core_members(mpt_high, detection=0, prevalence = 0.7)

# What are these ASVs? you can code it in two different ways to see the same things
tax_table(prune_taxa(low_ASVs,mpt_final))
tax_table(prune_taxa(int_ASVs,mpt_final))
tax_table(prune_taxa(high_ASVs,mpt_final))

# can plot those ASVs' relative abundance
prune_taxa(low_ASVs,mpt_RA) %>% 
  plot_bar(fill="Family") + 
  facet_wrap(.~`cn_category`, scales ="free")
prune_taxa(int_ASVs,mpt_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`cn_category`, scales ="free")
prune_taxa(high_ASVs,mpt_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`cn_category`, scales ="free")

# Notice that in this dataset, there are very few CORE microbiome members. This is common
### What if we wanted to make a Venn diagram of all the ASVs that showed up in each treatment?
low_list <- core_members(mpt_low, detection=0.001, prevalence = 0.10)
int_list <- core_members(mpt_int, detection=0.001, prevalence = 0.10)
high_list <- core_members(mpt_high, detection=0.001, prevalence = 0.10)

cn_list_full <- list(Low = low_list, Intermediate = int_list, High = high_list)

library("sf")

# Create a Venn diagram using all the ASVs shared and unique to antibiotic users and non users
first_venn <- ggVennDiagram(x = cn_list_full)

# Extract unique and overlapping ASVs
unique_low <- setdiff(cn_list_full[[1]], union(cn_list_full[[2]], cn_list_full[[3]]))
unique_int <- setdiff(cn_list_full[[2]], union(cn_list_full[[1]], cn_list_full[[3]]))
unique_high <- setdiff(cn_list_full[[3]], union(cn_list_full[[1]], cn_list_full[[2]]))

overlap_low_int <- intersect(cn_list_full[[1]], cn_list_full[[2]])
overlap_low_high <- intersect(cn_list_full[[1]], cn_list_full[[3]])
overlap_int_high <- intersect(cn_list_full[[2]], cn_list_full[[3]])
overlap_low_int_high <- Reduce(intersect, cn_list_full)

# Combine into a named list for easy iteration
asv_categories <- list(
  "unique_low" = unique_low,
  "unique_int" = unique_int,
  "unique_high" = unique_high,
  "overlap_low_int" = overlap_low_int,
  "overlap_low_high" = overlap_low_high,
  "overlap_int_high" = overlap_int_high,
  "overlap_low_int_high" = overlap_low_int_high
)

# Function to create a dataframe with ASV IDs and their taxonomic information
get_taxonomic_dataframe <- function(asv_list, phyloseq_obj) {
  if (length(asv_list) > 0) {
    # Prune the phyloseq object to only include the ASVs in the list
    pruned_phyloseq <- prune_taxa(asv_list, phyloseq_obj)
    # Extract taxonomic table and convert it to a dataframe
    tax_df <- as.data.frame(tax_table(pruned_phyloseq))
    # Add a column for ASV IDs
    tax_df <- tax_df %>% 
      mutate(ASV = rownames(tax_df)) %>% 
      select(ASV, everything())  # Make ASV the first column
    return(tax_df)
  } else {
    return(data.frame())  # Return an empty dataframe if no ASVs are present
  }
}

# Create a list to store dataframes for each category
taxonomic_dataframes <- list()

# Loop through each category and create the corresponding dataframe
for (category in names(asv_categories)) {
  cat("\nCreating dataframe for:", category, "\n")
  asvs <- asv_categories[[category]]
  taxonomic_dataframes[[category]] <- get_taxonomic_dataframe(asvs, mpt_final)
}

# Accessing each dataframe individually by name
unique_low_df <- taxonomic_dataframes[["unique_low"]]
unique_int_df <- taxonomic_dataframes[["unique_int"]]
unique_high_df <- taxonomic_dataframes[["unique_high"]]
overlap_low_int_df <- taxonomic_dataframes[["overlap_low_int"]]
overlap_low_high_df <- taxonomic_dataframes[["overlap_low_high"]]
overlap_int_high_df <- taxonomic_dataframes[["overlap_int_high"]]
overlap_low_int_high_df <- taxonomic_dataframes[["overlap_low_int_high"]]

# ISA #

library(tidyverse)
library(phyloseq)
library(indicspecies)

#### Indicator Species/Taxa Analysis ####
# glom to Genus
mpt_genus <- tax_glom(mpt_final, "Genus", NArm = FALSE)
mpt_genus_RA <- transform_sample_counts(mpt_genus, fun=function(x) x/sum(x))

#ISA
isa_mpt <- multipatt(t(otu_table(mpt_genus_RA)), cluster = sample_data(mpt_genus_RA)$`cn_category`)
summary(isa_mpt)
taxtable <- tax_table(mpt_final) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# consider that your table is only going to be resolved up to the genus level, be wary of 
# anything beyond the glomed taxa level
isa_df <- isa_mpt$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) 

# Extract ASVs for each cn_category based on the values in s.High, s.Intermediate, and s.Low
low_taxa <- isa_df %>% filter(s.Low == 1) %>% pull(ASV)
int_taxa <- isa_df %>% filter(s.Intermediate == 1) %>% pull(ASV)
high_taxa <- isa_df %>% filter(s.High == 1) %>% pull(ASV)

taxa_cn_full <- list(Low = low_taxa, Intermediate = int_taxa, High = high_taxa)

library("sf")

# Create a Venn diagram using all the ASVs shared and unique to categories in taxa_cn_full
second_venn <- ggVennDiagram(x = taxa_cn_full) +
  theme(
    text = element_text(size = 15)
  )
second_venn

ggsave("wetlands_venn_isa.png"
       , second_venn
       , height=8, width =8)

# Extract unique and overlapping ASVs
unique_low_second <- setdiff(taxa_cn_full[[1]], union(taxa_cn_full[[2]], taxa_cn_full[[3]]))
unique_int_second <- setdiff(taxa_cn_full[[2]], union(taxa_cn_full[[1]], taxa_cn_full[[3]]))
unique_high_second <- setdiff(taxa_cn_full[[3]], union(taxa_cn_full[[1]], taxa_cn_full[[2]]))

overlap_low_int_second <- intersect(taxa_cn_full[[1]], taxa_cn_full[[2]])
overlap_low_high_second <- intersect(taxa_cn_full[[1]], taxa_cn_full[[3]])
overlap_int_high_second <- intersect(taxa_cn_full[[2]], taxa_cn_full[[3]])
overlap_low_int_high_second <- Reduce(intersect, taxa_cn_full)

# Combine into a named list for easy iteration
asv_categories_second <- list(
  "unique_low_second" = unique_low_second,
  "unique_int_second" = unique_int_second,
  "unique_high_second" = unique_high_second,
  "overlap_low_int_second" = overlap_low_int_second,
  "overlap_low_high_second" = overlap_low_high_second,
  "overlap_int_high_second" = overlap_int_high_second,
  "overlap_low_int_high_second" = overlap_low_int_high_second
)

# Function to create a dataframe with ASV IDs and their taxonomic information
get_taxonomic_dataframe_second <- function(asv_list, phyloseq_obj) {
  if (length(asv_list) > 0) {
    # Prune the phyloseq object to only include the ASVs in the list
    pruned_phyloseq <- prune_taxa(asv_list, phyloseq_obj)
    # Extract taxonomic table and convert it to a dataframe
    tax_df <- as.data.frame(tax_table(pruned_phyloseq))
    # Add a column for ASV IDs
    tax_df <- tax_df %>% 
      mutate(ASV = rownames(tax_df)) %>% 
      select(ASV, everything())  # Make ASV the first column
    return(tax_df)
  } else {
    return(data.frame())  # Return an empty dataframe if no ASVs are present
  }
}

# Create a list to store dataframes for each category
taxonomic_dataframes_second <- list()

# Loop through each category and create the corresponding dataframe
for (category in names(asv_categories_second)) {
  cat("\nCreating dataframe for:", category, "\n")
  asvs <- asv_categories_second[[category]]
  taxonomic_dataframes_second[[category]] <- get_taxonomic_dataframe_second(asvs, mpt_final)
}

# Accessing each dataframe individually by name
unique_low_second_df <- taxonomic_dataframes_second[["unique_low_second"]]
unique_int_second_df <- taxonomic_dataframes_second[["unique_int_second"]]
if (nrow(unique_int_second_df) == 0) {
  # Taxonomic headers extracted from tax_table
  taxonomic_headers <- c("ASV", colnames(tax_table(mpt_final)))
  unique_int_second_df <- as.data.frame(matrix(ncol = length(taxonomic_headers), nrow = 0))
  colnames(unique_int_second_df) <- taxonomic_headers
}
unique_int_second_df <- unique_int_second_df %>%
  mutate(ASV = as.character(ASV))
unique_high_second_df <- taxonomic_dataframes_second[["unique_high_second"]]
overlap_low_int_second_df <- taxonomic_dataframes_second[["overlap_low_int_second"]]
if (nrow(overlap_low_int_second_df) == 0) {
  # Create a properly structured empty dataframe for overlap_low_int_second_df
  taxonomic_headers <- c("ASV", colnames(tax_table(mpt_final)))
  overlap_low_int_second_df <- as.data.frame(matrix(ncol = length(taxonomic_headers), nrow = 0))
  colnames(overlap_low_int_second_df) <- taxonomic_headers
}
overlap_low_int_second_df <- overlap_low_int_second_df %>%
  mutate(ASV = as.character(ASV))
overlap_low_high_second_df <- taxonomic_dataframes_second[["overlap_low_high_second"]]
if (nrow(overlap_low_high_second_df) == 0) {
  # Create a properly structured empty dataframe for overlap_low_high_second_df
  taxonomic_headers <- c("ASV", colnames(tax_table(mpt_final)))
  overlap_low_high_second_df <- as.data.frame(matrix(ncol = length(taxonomic_headers), nrow = 0))
  colnames(overlap_low_high_second_df) <- taxonomic_headers
}
overlap_low_high_second_df <- overlap_low_high_second_df %>%
  mutate(ASV = as.character(ASV))
overlap_int_high_second_df <- taxonomic_dataframes_second[["overlap_int_high_second"]]
overlap_low_int_high_second_df <- taxonomic_dataframes_second[["overlap_low_int_high_second"]]
if (nrow(overlap_low_int_high_second_df) == 0) {
  # Create a properly structured empty dataframe for overlap_low_int_high_second_df
  taxonomic_headers <- c("ASV", colnames(tax_table(mpt_final)))
  overlap_low_int_high_second_df <- as.data.frame(matrix(ncol = length(taxonomic_headers), nrow = 0))
  colnames(overlap_low_int_high_second_df) <- taxonomic_headers
}
overlap_low_int_high_second_df <- overlap_low_int_high_second_df %>%
  mutate(ASV = as.character(ASV))

# Perform an inner join on the ASV column to keep only common ASVs
unique_high <- inner_join(unique_high_df, unique_high_second_df, by = "ASV")
unique_int <- inner_join(unique_int_df, unique_int_second_df, by = "ASV")
unique_low <- inner_join(unique_low_df, unique_low_second_df, by = "ASV")
overlap_low_int <- inner_join(overlap_low_int_df, overlap_low_int_second_df, by = "ASV")
overlap_low_high <- inner_join(overlap_low_high_df, overlap_low_high_second_df, by = "ASV")
overlap_int_high <- inner_join(overlap_int_high_df, overlap_int_high_second_df, by = "ASV")
overlap_low_int_high <- inner_join(overlap_low_int_high_df, overlap_low_int_high_second_df, by = "ASV")

cn_full <- list(Low = unique_low$ASV, Intermediate = unique_int$ASV, High = unique_high$ASV)

library("sf")

# Create a Venn diagram using all the ASVs shared and unique to categories in taxa_cn_full
venn <- ggVennDiagram(x = cn_full)

write.csv(unique_low, "wetlands_taxa_low.csv", row.names = TRUE )

ggsave("wetlands_venn_combined.png"
       , venn
       , height=8, width =8)

ggsave("wetlands_venn_core.png"
       , first_venn
       , height=8, width =8)


