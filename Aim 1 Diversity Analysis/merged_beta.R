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
bc_dm <- phyloseq::distance(mpt_rare, method="bray")
uu_dm <- phyloseq::distance(mpt_rare, method="unifrac")
wu_dm <- phyloseq::distance(mpt_rare, method="wunifrac")

# create pcoa coordinate systen
pcoa_bc <- ordinate(mpt_rare, method="PCoA", distance=bc_dm)
pcoa_uu <- ordinate(mpt_rare, method="PCoA", distance=uu_dm)
pcoa_wu <- ordinate(mpt_rare, method="PCoA", distance=wu_dm)

# Set row names of metadata to sample IDs for cross-referencing
rownames(meta) <- meta$`#SampleID`

# Ensure metadata matches the distance matrices using row names
meta_filtered_bc <- meta[rownames(as.matrix(bc_dm)), , drop = FALSE]
meta_filtered_uu <- meta[rownames(as.matrix(uu_dm)), , drop = FALSE]
meta_filtered_wu <- meta[rownames(as.matrix(wu_dm)), , drop = FALSE]

# Run PERMANOVA for Bray-Curtis distance matrix using adonis2
permanova_bc <- adonis2(bc_dm ~ cn_category_per_dataset, data = meta_filtered_bc, permutations = 999)
print(permanova_bc)

# Run PERMANOVA for Unweighted UniFrac distance matrix using adonis2
permanova_uu <- adonis2(uu_dm ~ cn_category_per_dataset, data = meta_filtered_uu, permutations = 999)
print(permanova_uu)

# Run PERMANOVA for Weighted UniFrac distance matrix using adonis2
permanova_wu <- adonis2(wu_dm ~ cn_category_per_dataset, data = meta_filtered_wu, permutations = 999)
print(permanova_wu)

# plot coordination
gg_pcoa_bc <- plot_ordination(mpt_rare, pcoa_bc, color = "cn_category_per_dataset", shape="subject") +
  labs(pch="Subject #", col = "C:N Category")
gg_pcoa_bc <- gg_pcoa_bc +
  stat_ellipse(aes(group = cn_category_per_dataset), level = 0.95, linetype = "dashed") +
  annotate("text", x = 0.3, y = -0.375, label = "p < 0.001 ***", size = 5)
gg_pcoa_bc

gg_pcoa_uu <- plot_ordination(mpt_rare, pcoa_uu, color = "cn_category_per_dataset", shape="subject") +
  labs(pch="Subject #", col = "C:N Category")
gg_pcoa_uu <- gg_pcoa_uu +
  stat_ellipse(aes(group = cn_category_per_dataset), level = 0.95, linetype = "dashed") +
  annotate("text", x = 0.4, y = -0.375, label = "p < 0.001 ***", size = 5)
gg_pcoa_uu

gg_pcoa_wu <- plot_ordination(mpt_rare, pcoa_wu, color = "cn_category_per_dataset", shape="subject") +
  labs(pch="Subject #", col = "C:N Category")
gg_pcoa_wu <- gg_pcoa_wu +
  stat_ellipse(aes(group = cn_category_per_dataset), level = 0.95, linetype = "dashed") +
  annotate("text", x = 0.0872, y = -0.04, label = "p < 0.001 ***", size = 5)
gg_pcoa_wu

# create box plots 
library(reshape2)
library(dplyr) 

# Bray-Curtis Distance Matrix
bray_dist_long <- melt(as.matrix(bc_dm))
colnames(bray_dist_long) <- c("Sample1", "Sample2", "Distance")
# Merge with meta to add 'cn_category_per_dataset' grouping information
bray_dist_long <- bray_dist_long %>%
  left_join(meta, by = c("Sample1" = "#SampleID")) %>%
  rename(Group = cn_category_per_dataset)

# Unweighted UniFrac Distance Matrix
unweighted_unifrac_long <- melt(as.matrix(uu_dm))
colnames(unweighted_unifrac_long) <- c("Sample1", "Sample2", "Distance")
unweighted_unifrac_long <- unweighted_unifrac_long %>%
  left_join(meta, by = c("Sample1" = "#SampleID")) %>%
  rename(Group = cn_category_per_dataset)

# Weighted UniFrac Distance Matrix
weighted_unifrac_long <- melt(as.matrix(wu_dm))
colnames(weighted_unifrac_long) <- c("Sample1", "Sample2", "Distance")
weighted_unifrac_long <- weighted_unifrac_long %>%
  left_join(meta, by = c("Sample1" = "#SampleID")) %>%
  rename(Group = cn_category_per_dataset)

library(ggplot2)

# Box Plot for Bray-Curtis Distance Matrix
gg_box_bc <- ggplot(bray_dist_long, aes(x = Group, y = Distance, fill = Group)) +
  geom_boxplot() +
  labs(title = "Beta Diversity Box Plot (Bray-Curtis) by cn_category_per_dataset",
       x = "Sample Group (cn_category_per_dataset)",
       y = "Beta Diversity (Bray-Curtis)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# Box Plot for Unweighted UniFrac Distance Matrix
gg_box_uu <- ggplot(unweighted_unifrac_long, aes(x = Group, y = Distance, fill = Group)) +
  geom_boxplot() +
  labs(title = "Beta Diversity Box Plot (Unweighted UniFrac) by cn_category_per_dataset",
       x = "Sample Group (cn_category_per_dataset)",
       y = "Beta Diversity (Unweighted UniFrac)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Box Plot for Weighted UniFrac Distance Matrix
gg_box_wu <- ggplot(weighted_unifrac_long, aes(x = Group, y = Distance, fill = Group)) +
  geom_boxplot() +
  labs(title = "Beta Diversity Box Plot (Weighted UniFrac) by cn_category_per_dataset",
       x = "Sample Group (cn_category_per_dataset)",
       y = "Beta Diversity (Weighted UniFrac)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plots
print(gg_box_bc)
print(gg_box_uu)
print(gg_box_wu)

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
  labs(title = "Beta Diversity by C:N category across Datasets",
       x = "C:N category",
       y = "Beta Diversity (Weighted UniFrac)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 20, face = "bold"),  
    axis.text = element_text(size = 14),                  
    axis.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),                
    legend.title = element_text(size = 18, face = "bold")
  ) +
  geom_signif(comparisons = list(
    c("High-Soil", "High-Wetlands"), 
    c("Intermediate-Soil", "Intermediate-Wetlands"), 
    c("Low-Soil", "Low-Wetlands")),
    annotations = "****",
    y_position = c(0.125, 0.125, 0.125),  # Setting the height for the bars
    tip_length = 0.01,
    textsize = 5)

# Print the plot
print(gg_box_wu_s)

ggsave("merged_wu.png"
       , gg_box_wu_s
       , height=8, width =12)

# Extract the fill colors from the merged plot
merged_colors <- ggplot_build(gg_box_wu_s)$data[[1]]$fill
unique_colors <- unique(merged_colors)
print(unique_colors)

# DESeq #
library(DESeq2)

mpt_plus1 <- transform_sample_counts(mpt_final, function(x) x+1)
mpt_deseq <- phyloseq_to_deseq2(mpt_plus1, ~`cn_category_per_dataset`)
DESEQ_mpt <- DESeq(mpt_deseq)

## Low ##
res_low <- results(DESEQ_mpt, tidy=TRUE, 
                   contrast = c("cn_category_per_dataset","Low-Wetlands","Low-Soil"))

vol_plot_low <- res_low %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
vol_plot_low

sigASVs_low <- res_low %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

sigASVs_vec_low <- sigASVs_low %>%
  pull(ASV)

mpt_DESeq_low <- prune_taxa(sigASVs_vec_low,mpt_final)
sigASVs_low <- tax_table(mpt_DESeq_low) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_low) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

log2_low <- ggplot(sigASVs_low) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
log2_low

# Count the number of positive log2FoldChange values
positive_count_low <- sum(sigASVs_low$log2FoldChange > 0)

# Count the number of negative log2FoldChange values
negative_count_low <- sum(sigASVs_low$log2FoldChange < 0)

# Create a new data frame 
count_low <- data.frame(
  log2FoldChange = c("Positive", "Negative"),
  Count = c(positive_count_low, negative_count_low)
)

## Int ##
res_int <- results(DESEQ_mpt, tidy=TRUE, 
                   contrast = c("cn_category_per_dataset","Intermediate-Wetlands","Intermediate-Soil"))

vol_plot_int <- res_int %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
vol_plot_int

sigASVs_int <- res_int %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

sigASVs_vec_int <- sigASVs_int %>%
  pull(ASV)

mpt_DESeq_int <- prune_taxa(sigASVs_vec_int,mpt_final)
sigASVs_int <- tax_table(mpt_DESeq_int) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_int) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

log2_int <- ggplot(sigASVs_int) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
log2_int

# Count the number of positive log2FoldChange values
positive_count_int <- sum(sigASVs_int$log2FoldChange > 0)

# Count the number of negative log2FoldChange values
negative_count_int <- sum(sigASVs_int$log2FoldChange < 0)

# Create a new data frame 
count_int <- data.frame(
  log2FoldChange = c("Positive", "Negative"),
  Count = c(positive_count_int, negative_count_int)
)

## High
res_high <- results(DESEQ_mpt, tidy=TRUE, 
                   contrast = c("cn_category_per_dataset","High-Wetlands","High-Soil"))

vol_plot_high <- res_high %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
vol_plot_high

sigASVs_high <- res_high %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

sigASVs_vec_high <- sigASVs_high %>%
  pull(ASV)

mpt_DESeq_high <- prune_taxa(sigASVs_vec_high,mpt_final)
sigASVs_high <- tax_table(mpt_DESeq_high) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_high) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

log2_high <- ggplot(sigASVs_high) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
log2_high

# Count the number of positive log2FoldChange values
positive_count_high <- sum(sigASVs_high$log2FoldChange > 0)

# Count the number of negative log2FoldChange values
negative_count_high <- sum(sigASVs_high$log2FoldChange < 0)

# Create a new data frame 
count_high <- data.frame(
  log2FoldChange = c("Positive", "Negative"),
  Count = c(positive_count_high, negative_count_high)
)



