# Load in Libraries
library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library(picante)
library(ggsignif)

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


#### Accessor functions ####
# These functions allow you to see or summarise data

# If we look at sample variables and decide we only want some variables, we can view them like so:
sample_variables(mpt)
# colnames(sample_data(atamaca))
get_variable(mpt, c("cn_category")) # equivalent to "select" in tidyverse

## Let's say we want to filter OTU table by sample. 
# We can first view sample names:
sample_names(mpt)
# How many samples do we have?
nsamples(mpt)
# What is the sum of reads in each sample?
sample_sums(mpt)
# Save the sample names of the 3 samples with the most reads
getsamps <- names(sort(sample_sums(mpt), decreasing = TRUE)[1:3])
# filter to see taxa abundances for each sample
get_taxa(mpt, getsamps) 

## Conversely, let's say we want to compare OTU abundance of most abundant OTU across samples
# Look at taxa names
taxa_names(mpt)
# How many taxa do we have? (ASVs)
ntaxa(mpt)
# What is the total read count for each taxa?
taxa_sums(mpt)
# Let's find the top 3 most abundant taxa
gettaxa <- names(sort(taxa_sums(mpt), decreasing = TRUE)[1:3] )
get_sample(mpt, gettaxa)


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
rarecurve(t(as.data.frame(otu_table(mpt_final))), cex=0.1)

## CODE RUN TILL ABOVE AND VERIFIED - WORK ON BELOW 
## Wetlands has not been run or verified at all yet - need to do still

mpt_rare <- rarefy_even_depth(mpt_final, rngseed = 8, sample.size = 2500)

#### Alpha diversity ######
plot_richness(mpt_rare) 

plot_richness(mpt_rare, measures = c("Shannon","Chao1")) 

gg_richness <- plot_richness(mpt_rare, x = "cn_category", measures = c("Shannon", "Observed", "Fisher")) +
  xlab("C:N Category") +
  geom_boxplot()
gg_richness

ggsave(filename = "plot_richness.png"
       , gg_richness
       , height=4, width=6)

estimate_richness(mpt_rare)

# phylogenetic diversity

# calculate Faith's phylogenetic diversity as PD
phylo_dist <- pd(t(otu_table(mpt_rare)), phy_tree(mpt_rare),
                 include.root=F) 

# add PD to metadata table
sample_data(mpt_rare)$PD <- phylo_dist$PD

# plot any metadata category against the PD
plot.pd <- ggplot(sample_data(mpt_rare), aes(cn_category, PD)) + 
  geom_boxplot() +
  xlab("C:N Category") +
  ylab("Phylogenetic Diversity")

# view plot
plot.pd

# Need to extract information
alphadiv <- estimate_richness(mpt_rare)
samp_dat <- sample_data(mpt_rare)
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)

View(samp_dat_wdiv)

allCounts <- as.vector(otu_table(mpt_rare))
allCounts <- allCounts[allCounts>0]
hist(allCounts)
hist(log(allCounts))

ggplot(samp_dat_wdiv, aes(x=`cn_category`, y=Shannon)) +
  geom_boxplot() +
  geom_point()

kruskal.test(Shannon ~ `cn_category`, data = samp_dat_wdiv)
kruskal.test(Observed ~ `cn_category`, data = samp_dat_wdiv)
kruskal.test(Fisher ~ `cn_category`, data = samp_dat_wdiv)

lm_ob_vs_site_log <- lm(log(Shannon) ~ `cn_category`, data=samp_dat_wdiv)
anova_ob_vs_site_log <- aov(lm_ob_vs_site_log)
summary(anova_ob_vs_site_log)
TukeyHSD(anova_ob_vs_site_log)

# Mapping the significance to a ggplot
Shannon <- ggplot(samp_dat_wdiv, aes(x=`cn_category`, y=Shannon)) +
  geom_boxplot() +
  geom_signif(comparisons = list(c("Low","High"), c("Low", "Intermediate")),
              y_position = c(6, 6.5),
              annotations = c("0.0002","0.0007"))
Shannon

lm_ob_vs_site_log_pd <- lm(log(PD) ~ `cn_category`, data=samp_dat_wdiv)
anova_ob_vs_site_log_pd <- aov(lm_ob_vs_site_log_pd)
summary(anova_ob_vs_site_log_pd)
TukeyHSD(anova_ob_vs_site_log_pd)

PD <- Shannon <- ggplot(samp_dat_wdiv, aes(x=`cn_category`, y=PD)) +
  geom_boxplot() +
  geom_signif(comparisons = list(c("Low","High"), c("Low", "Intermediate"), c("Intermediate", "High")),
              y_position = c(70, 65, 60),
              annotations = c("0.0000000","0.004", "0.0007"))
PD

