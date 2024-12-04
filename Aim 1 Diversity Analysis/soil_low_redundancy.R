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
# Filter meta DataFrame to include ONLY "Low in cn_category
meta <- meta[meta$cn_category %in% c("Low"), ]
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

# More libraries #
library(readr)
library(data.table)

#### Import files and preparing tables ####
# Importing the pathway PICrsut2
abundance_file <- "../qiime/Soil/picrust2/pathway_abundance.tsv"
abundance_data <- fread(abundance_file, sep = "\t", header = TRUE, strip.white = TRUE)
abundance_data  =as.data.frame(abundance_data)

setnames(abundance_data, old = "#OTU ID", new = "pathway")

pathway_abundance <- abundance_data

# Set pathway names as rownames and remove the first column
rownames(pathway_abundance) <- pathway_abundance$pathway
pathway_abundance <- pathway_abundance[, -1]
pathway_abundance <- t(pathway_abundance)

# Match samples between phyloseq and pathway_abundance
sample_ids_physeq <- sample_names(mpt)
pathway_abundance <- pathway_abundance[rownames(pathway_abundance) %in% sample_ids_physeq, ]
pathway_abundance <- pathway_abundance[sample_ids_physeq, ]

# Calculate functional distances using Bray-Curtis distance
dis <- vegdist(pathway_abundance, method = "bray")

# Normalize the distances to ensure they are between 0 and 1
dis <- dis / max(dis)

# Assuming your phyloseq object is called 'physeq'
otu_table_data <- otu_table(mpt)
# Convert the OTU table to a data frame
comm <- as.data.frame(otu_table_data)
comm <- t(comm)

dis <- vegdist(comm, method = "bray")

# Normalize the distances to ensure they are between 0 and 1
dis <- dis / max(dis)

## INPUT uniquenesss code ##
uniqueness <- function(comm, dis, tol = 1e-8, abundance = TRUE){
  
  if(!is.null(colnames(comm)) & !is.null(attributes(dis)$Labels)){
    if(any(!colnames(comm)%in%attributes(dis)$Labels)) stop("One or several species in comm are not in dis; check species names in comm and in dis")
    else dis <- as.dist(as.matrix(dis)[colnames(comm), colnames(comm)])
  }
  else if(ncol(comm)!=attributes(dis)$Size) stop("the number of species in comm must be equal to that in dis") 		
  
  D <- as.matrix(dis)
  if(!abundance){
    comm[comm>0] <- 1
  }
  commt <- as.data.frame(t(comm))
  
  funK <- function(v){
    p <- v/sum(v)
    K <- apply(D, 1, function(x) sum(x*p))
    K[p<tol] <- NA
    return(K)
  }
  V <- cbind.data.frame(sapply(commt, funK))
  rownames(V) <- colnames(comm)
  colnames(V) <- rownames(comm)
  funKbar <- function(v){
    p <- v/sum(v)
    Kbar <- sapply(1:nrow(D), function(i) sum(D[i,]*v/sum(v[-i])))
    Kbar[p<tol] <- NA
    return(Kbar)
  }
  Kbar <- cbind.data.frame(sapply(commt, funKbar))
  rownames(Kbar) <- colnames(comm)
  colnames(Kbar) <- rownames(comm)
  funQ <- function(v){
    p <- v/sum(v)
    Q <- t(p)%*%D%*%p
    return(Q)
  }
  Q <- unlist(sapply(commt, funQ))
  
  funSim <- function(v){
    p <- v/sum(v)
    S <- 1-sum(p^2)
    return(S)
  }
  Sim <- unlist(sapply(commt, funSim))
  
  funN <- function(v){
    p <- v/sum(v)
    N <- length(p[p>tol])
    return(N)
  }
  N <- unlist(sapply(commt, funN))
  U <- Q/Sim
  
  
  red <- cbind.data.frame(N=N, D=Sim, Q=Q, U=U)
  rownames(red) <- rownames(comm)
  
  res <- list()
  res$Kbar <- Kbar
  res$V <- V
  res$red <- red
  return(res)
  
}

result <- uniqueness(comm, dis, tol = 1e-8, abundance = TRUE)

