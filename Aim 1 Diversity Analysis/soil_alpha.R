# Load in Libraries
install.packages("ggthemes")
library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library(picante)
library(ggsignif)
library(ggthemes)

#### Load data ####
# Change file paths as necessary
metafp <- "../wrangling/soil_metadata.tsv"
meta <- read_delim(metafp, delim="\t")

otufp <- "../qiime/Soil/soil_export/table_export/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "../qiime/Soil/soil_export/taxonomy_export/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "../qiime/Soil/soil_export/rooted_tree_export/tree.nwk"
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
# Filter meta DataFrame to exclude "Very Low" and "Very High" in cn_category
meta <- meta[!meta$cn_category %in% c("Very Low", "Very High"), ]
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
mpt_final <- subset_samples(mpt_filt_nolow_samps, !is.na(cn_category) )

# Rarefy samples
# rngseed sets a random number. If you want to reproduce this exact analysis, you need
# to set rngseed the same number each time
# t transposes the table to use rarecurve function
# cex decreases font size
# rarecurve(t(as.data.frame(otu_table(mpt_final))), cex=0.1)

mpt_rare <- rarefy_even_depth(mpt_final, rngseed = 8, sample.size = 2500)
plot_richness(mpt_rare)

save(mpt_final, file = "soil_final.RData")
save(mpt_rare, file = "soil_rare.RData")

#### Alpha diversity ######
# phylogenetic diversity

# calculate Faith's phylogenetic diversity as PD
phylo_dist <- pd(t(otu_table(mpt_rare)), phy_tree(mpt_rare),
                 include.root=F) 

# add PD to metadata table
sample_data(mpt_rare)$PD <- phylo_dist$PD

# Need to extract information
alphadiv <- estimate_richness(mpt_rare)
samp_dat <- sample_data(mpt_rare)
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)

kruskal.test(Shannon ~ `cn_category`, data = samp_dat_wdiv)
kruskal.test(Observed ~ `cn_category`, data = samp_dat_wdiv)
kruskal.test(Fisher ~ `cn_category`, data = samp_dat_wdiv)
kruskal.test(PD ~ `cn_category`, data = samp_dat_wdiv)

lm_ob_vs_site_log_pd <- lm(log(PD) ~ `cn_category`, data=samp_dat_wdiv)
anova_ob_vs_site_log_pd <- aov(lm_ob_vs_site_log_pd)
summary(anova_ob_vs_site_log_pd)
TukeyHSD(anova_ob_vs_site_log_pd)

samp_dat_wdiv$cn_category <- factor(samp_dat_wdiv$cn_category, levels = c("Low", "Intermediate", "High"))

PD <- ggplot(samp_dat_wdiv, aes(x=`cn_category`, y=PD, fill = cn_category)) +
  geom_boxplot() +
  labs(x="C:N Category", y="Faith's PD", fill = expression(bold("C:N Category"))) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 25, face = "bold"),  
    axis.text = element_text(size = 25),                  
    axis.title = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 20),                
    legend.title = element_text(size = 25, face = "bold")) +
  geom_signif(comparisons = list(c("Low","High"), c("Low", "Intermediate"), c("Intermediate", "High")),
              y_position = c(70, 65, 60),
              annotations = c("****","**", "***"),
              textsize = 8) + 
  scale_fill_manual(values = c("Low" = "#E69F00", 
                               "Intermediate" = "#56B4E9", 
                               "High" = "#009E73"))

PD

samp_dat_wdiv$cn_category <- factor(samp_dat_wdiv$cn_category, levels = c("Low", "Intermediate", "High"))

lm_ob_vs_site_log_shannon <- lm(log(Shannon) ~ `cn_category`, data=samp_dat_wdiv)
anova_ob_vs_site_log_shannon <- aov(lm_ob_vs_site_log_shannon)
summary(anova_ob_vs_site_log_shannon)
TukeyHSD(anova_ob_vs_site_log_shannon)

Shannon <- ggplot(samp_dat_wdiv, aes(x=`cn_category`, y=Shannon, fill = cn_category)) +
  geom_boxplot() +
  labs(x="C:N Category", y="Shannon Evenness", fill = expression(bold("C:N Category"))) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 25, face = "bold"),  
    axis.text = element_text(size = 25),                  
    axis.title = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 20),                
    legend.title = element_text(size = 25, face = "bold")
  ) +
  geom_signif(comparisons = list(c("Low","High"), c("Low", "Intermediate")),
              y_position = c(6.25, 5.9),
              annotations = c("***","***"),
              textsize = 8) + 
  scale_fill_manual(values = c("Low" = "#E69F00", 
                               "Intermediate" = "#56B4E9", 
                               "High" = "#009E73"))
Shannon

ggsave("soil_PD.png"
       , PD
       , height=8, width =12)

ggsave("soil_shannon.png"
       , Shannon
       , height=8, width =12)

PD2 <- ggplot(samp_dat_wdiv, aes(x=`cn_category`, y=PD, fill = cn_category)) +
  geom_boxplot() +
  labs(x="C:N Category", y="Faith's PD", fill = expression(bold("C:N Category"))) +
  theme(
    plot.title = element_text(size = 35, face = "bold"),  
    axis.text = element_text(size = 25),                  
    axis.title = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 20),                
    legend.title = element_text(size = 25, face = "bold")
  ) +
  geom_signif(comparisons = list(c("Low","High"), c("Low", "Intermediate"), c("Intermediate", "High")),
              y_position = c(70, 65, 60),
              annotations = c("****","**", "***"),
              textsize = 8) + 
  scale_fill_manual(values = c("Low" = "#619CFF", 
                               "Intermediate" = "#00BA38", 
                               "High" = "#F8766D"))
PD2

Shannon2 <- ggplot(samp_dat_wdiv, aes(x=`cn_category`, y=Shannon, fill = cn_category)) +
  geom_boxplot() +
  labs(x="C:N Category", y="Shannon Evenness", title="Forest Soil", fill = expression(bold("C:N Category"))) +
  theme(
    plot.title = element_text(size = 35, face = "bold"),  
    axis.text = element_text(size = 25),                  
    axis.title = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 20),                
    legend.title = element_text(size = 25, face = "bold")
  ) +
  geom_signif(comparisons = list(c("Low","High"), c("Low", "Intermediate")),
              y_position = c(6.25, 5.9),
              annotations = c("***","***"),
              textsize = 8) + 
  scale_fill_manual(values = c("Low" = "#619CFF", 
                               "Intermediate" = "#00BA38", 
                               "High" = "#F8766D"))
Shannon2

ggsave("soil_PD2.png"
       , PD2
       , height=8, width =12)
ggsave("soil_shannon2.png"
       , Shannon2
       , height=8, width =12)

# DESeq #
library(DESeq2)

mpt_plus1 <- transform_sample_counts(mpt_final, function(x) x+1)
mpt_deseq <- phyloseq_to_deseq2(mpt_plus1, ~`cn_category`)
DESEQ_mpt <- DESeq(mpt_deseq)

## Low vs High ##
res_lvh <- results(DESEQ_mpt, tidy=TRUE, 
               contrast = c("cn_category","Low","High"))

vol_plot_lvh <- res_lvh %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
vol_plot_lvh

sigASVs_lvh <- res_lvh %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

sigASVs_vec_lvh <- sigASVs_lvh %>%
  pull(ASV)

mpt_DESeq_lvh <- prune_taxa(sigASVs_vec_lvh,mpt_final)
sigASVs_lvh <- tax_table(mpt_DESeq_lvh) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_lvh) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

log2_lvh <- ggplot(sigASVs_lvh) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
log2_lvh

# Count the number of positive log2FoldChange values
positive_count_lvh <- sum(sigASVs_lvh$log2FoldChange > 0)

# Count the number of negative log2FoldChange values
negative_count_lvh <- sum(sigASVs_lvh$log2FoldChange < 0)

# Create a new data frame 
count_lvh <- data.frame(
  log2FoldChange = c("Positive", "Negative"),
  Count = c(positive_count_lvh, negative_count_lvh)
)

## Int vs High ##
res_ivh <- results(DESEQ_mpt, tidy=TRUE, 
                   contrast = c("cn_category","Intermediate","High"))

vol_plot_ivh <- res_ivh %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
vol_plot_ivh

sigASVs_ivh <- res_ivh %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

sigASVs_vec_ivh <- sigASVs_ivh %>%
  pull(ASV)

mpt_DESeq_ivh <- prune_taxa(sigASVs_vec_ivh,mpt_final)
sigASVs_ivh <- tax_table(mpt_DESeq_ivh) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_ivh) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

log2_ivh <- ggplot(sigASVs_ivh) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
log2_ivh

# Count the number of positive log2FoldChange values
positive_count_ivh <- sum(sigASVs_ivh$log2FoldChange > 0)

# Count the number of negative log2FoldChange values
negative_count_ivh <- sum(sigASVs_ivh$log2FoldChange < 0)

# Create a new data frame 
count_ivh <- data.frame(
  log2FoldChange = c("Positive", "Negative"),
  Count = c(positive_count_ivh, negative_count_ivh)
)

## Low vs Int
res_lvi <- results(DESEQ_mpt, tidy=TRUE, 
                   contrast = c("cn_category","Low","Intermediate"))

vol_plot_lvi <- res_lvi %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
vol_plot_lvi

sigASVs_lvi <- res_lvi %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

sigASVs_vec_lvi <- sigASVs_lvi %>%
  pull(ASV)

mpt_DESeq_lvi <- prune_taxa(sigASVs_vec_lvi,mpt_final)
sigASVs_lvi <- tax_table(mpt_DESeq_lvi) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_lvi) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

log2_lvi <- ggplot(sigASVs_lvi) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
log2_lvi

# Count the number of positive log2FoldChange values
positive_count_lvi <- sum(sigASVs_lvi$log2FoldChange > 0)

# Count the number of negative log2FoldChange values
negative_count_lvi <- sum(sigASVs_lvi$log2FoldChange < 0)

# Create a new data frame 
count_lvi <- data.frame(
  log2FoldChange = c("Positive", "Negative"),
  Count = c(positive_count_lvi, negative_count_lvi)
)


