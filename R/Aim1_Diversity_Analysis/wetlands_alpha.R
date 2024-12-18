# Load in Libraries
library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library(picante)
library(ggsignif)
library(ggthemes)
#### Load data ####
# Change file paths as necessary
metafp <- "../Wrangling/wetlands_metadata.tsv"
meta <- read_delim(metafp, delim="\t")

otufp <- "../../Qiime/Wetlands/wetlands_export/table_export/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "../../Qiime/Wetlands/wetlands_export/taxonomy_export/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "../../Qiime/Wetlands/wetlands_export/rooted_tree_export/tree.nwk"
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
mpt_rare <- rarefy_even_depth(mpt_final, rngseed = 8, sample.size = 20000)

save(mpt_final, file = "wetlands_final.RData")
save(mpt_rare, file = "wetlands_rare.RData")

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
  labs(x="C:N Category", y="Faith's PD",fill = expression(bold("C:N Category"))) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 25, face = "bold"),  
    axis.text = element_text(size = 25),                  
    axis.title = element_text(size = 25, face = "bold"),
    legend.position = "none") +
  geom_signif(comparisons = list(c("Low","High")),
              y_position = c(175),
              annotations = c("**"),
              textsize = 8) + 
  scale_fill_manual(values = c("Low" = "#E69F00", 
                               "Intermediate" = "#56B4E9", 
                               "High" = "#009E73")) +
  scale_y_continuous(limits = c(10, NA))

PD

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
    legend.position = "none") +
  scale_fill_manual(values = c("Low" = "#E69F00", 
                               "Intermediate" = "#56B4E9", 
                               "High" = "#009E73")) +
  scale_y_continuous(limits = c(3, NA))
Shannon

ggsave("wetlands_PD.png"
       , PD
       , height=8, width =12)

ggsave("wetlands_shannon.png"
       , Shannon
       , height=8, width =12)
