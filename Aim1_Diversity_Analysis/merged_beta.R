# Load in Libraries
library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library(picante)
library(ggsignif)

#### Load data ####
# Change file paths as necessary
metafp <- "../wrangling/merged_wrangled_metadata.tsv"
meta <- read_delim(metafp, delim="\t")

otufp <- "../qiime/Merged/merged_export/table_export/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "../qiime/Merged/merged_export/taxonomy_export/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "../qiime/Merged/merged_export/rooted_tree_export/tree.nwk"
phylotree <- read.tree(phylotreefp)

#### Format OTU table ####
# OTU tables should be a matrix
# with rownames and colnames as OTUs and sampleIDs, respectively
# Note: tibbles do not allow rownames so if you imported with read_delim, change back

# save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

#### Format sample metadata ####
# Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1])
# Make sampleids the rownames
rownames(samp_df)<- meta$'#SampleID'
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)

#### Formatting taxonomy ####
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)

#### Create phyloseq object ####
# Merge all into a phyloseq object
mpt <- phyloseq(OTU, SAMP, TAX, phylotree)

#### Looking at phyloseq object #####
# View components of phyloseq object with the following commands
otu_table(mpt)
sample_data(mpt)
tax_table(mpt)
phy_tree(mpt)

######### ANALYZE ##########
# Remove non-bacterial sequences, if any
mpt_filt <- subset_taxa(mpt,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# Remove ASVs that have less than 5 counts total
mpt_filt_nolow <- filter_taxa(mpt_filt, function(x) sum(x)>5, prune = TRUE)
# Remove samples with less than 100 reads - poor sequencing run
mpt_filt_nolow_samps <- prune_samples(sample_sums(mpt_filt_nolow)>100, mpt_filt_nolow)
# Remove samples where cn_category is na
mpt_final <- subset_samples(mpt_filt_nolow_samps, !is.na(cn_category_per_dataset) )

# Rarefy samples
# rngseed sets a random number. If you want to reproduce this exact analysis, you need
# to set rngseed the same number each time
# t transposes the table to use rarecurve function
# cex decreases font size
# rarecurve(t(as.data.frame(otu_table(mpt_final))), cex=0.1)

mpt_rare <- rarefy_even_depth(mpt_final, rngseed = 8, sample.size = 5000)

#### Beta diversity #####
wu_dm <- phyloseq::distance(mpt_rare, method="wunifrac")

# create pcoa coordinate systen
pcoa_wu <- ordinate(mpt_rare, method="PCoA", distance=wu_dm)

# Set row names of metadata to sample IDs for cross-referencing
rownames(meta) <- meta$`#SampleID`

# Ensure metadata matches the distance matrices using row names
meta_filtered_wu <- meta[rownames(as.matrix(wu_dm)), , drop = FALSE]

# Run PERMANOVA for Weighted UniFrac distance matrix using adonis2
permanova_wu <- adonis2(wu_dm ~ cn_category_per_dataset, data = meta_filtered_wu, permutations = 999)
print(permanova_wu)

# plot coordination
gg_pcoa_wu <- plot_ordination(mpt_rare, pcoa_wu, color = "cn_category_per_dataset", shape="subject") +
  labs(pch="Subject #", col = "C:N Category")
gg_pcoa_wu_s <- plot_ordination(mpt_rare, pcoa_wu, color = "cn_category_per_dataset", shape="subject") +
  labs(pch="Subject #", col = "C:N Category") +
  theme(
  plot.title = element_text(size = 25, face = "bold"),  
  axis.text = element_text(size = 25),                  
  axis.title = element_text(size = 25, face = "bold"),
  legend.text = element_text(size = 20),                
  legend.title = element_text(size = 25, face = "bold")) +
  stat_ellipse(aes(group = cn_category_per_dataset), level = 0.95, linetype = "dashed") +
  annotate("text", x = 0.0872, y = -0.04, label = "p < 0.001", size = 5, fontface = "bold") +
  annotate("text", x = 0.047, y = 0.035, label = "Wetlands Soil", size = 5, fontface = "bold") +
 annotate("text", x = -0.016, y = 0.003 , label = "Forest Soil", size = 5, fontface = "bold")
gg_pcoa_wu_s

ggsave("merged_pcoa.png"
       , gg_pcoa_wu_s
       , height=8, width =12)

# create box plots 
library(reshape2)
library(dplyr) 

# Weighted UniFrac Distance Matrix
weighted_unifrac_long <- melt(as.matrix(wu_dm))
colnames(weighted_unifrac_long) <- c("Sample1", "Sample2", "Distance")
weighted_unifrac_long <- weighted_unifrac_long %>%
  left_join(meta, by = c("Sample1" = "#SampleID")) %>%
  rename(Group = cn_category_per_dataset)

library(ggplot2)

# Run Kruskal-Wallis test on the Weighted UniFrac distances
kruskal_test_wu <- kruskal.test(Distance ~ Group, data = weighted_unifrac_long)
kruskal_test_wu

library(dunn.test)

# Run Dunn's post-hoc test after the Kruskal-Wallis test
dunn_test_wu <- dunn.test(weighted_unifrac_long$Distance, weighted_unifrac_long$Group, method = "bonferroni")

library(ggsignif)

# Box Plot for Weighted UniFrac Distance Matrix Grouped by cn_category_per_dataset
gg_box_wu_s <- ggplot(weighted_unifrac_long, aes(x = Group, y = Distance, fill = Group)) +
  geom_boxplot() +
  labs(x = "C:N category",
       y = "Beta Diversity (Weighted UniFrac)",
      fill = expression(bold("C:N Category"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  
        axis.text = element_text(size = 16),                  
        axis.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 20),                
        legend.title = element_text(size = 25, face = "bold")
  ) +
  geom_signif(comparisons = list(
    c("High-Soil", "High-Wetlands"), 
    c("Intermediate-Soil", "Intermediate-Wetlands"), 
    c("Low-Soil", "Low-Wetlands")),
    annotations = "****",
    y_position = c(0.125, 0.125, 0.125),  # Setting the height for the bars
    tip_length = 0.01,
    textsize = 5)
gg_box_wu_s

ggsave("merged_wu.png"
       , gg_box_wu_s
       , height=8, width =12)
