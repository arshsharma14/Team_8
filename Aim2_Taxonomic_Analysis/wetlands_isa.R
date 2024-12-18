# Load Libraries #
library(tidyverse)
library(phyloseq)
library(indicspecies)
library(ggVennDiagram)

#### Load data ####
load("../Aim1_Diversity_Analysis/wetlands_final.RData")

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

