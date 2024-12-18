# Load Libraries #
library(tidyverse)
library(phyloseq)
library(DESeq2)

#### Load data ####
load("../Aim1_Diversity_Analysis/soil_final.RData")

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

