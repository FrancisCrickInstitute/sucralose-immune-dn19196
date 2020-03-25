library(tidyverse)
library(dada2)
library(usethis)
library(openxlsx)

dat <- read.table("results/abundance_table/filtered/feature-table.tsv", header=TRUE, comment.char="", skip=1, sep="\t") %>%
  tibble::column_to_rownames(var="X.OTU.ID") %>%
  as.matrix() %>%
  otu_table(taxa_are_rows=TRUE)


tax <- read.table("results/taxonomy/taxonomy.tsv", header=TRUE, comment.char="", sep="\t") %>%
  dplyr::select(-Confidence) %>%
  dplyr::mutate(Taxon=gsub("D_[0-6]__","", Taxon)) %>%
  tibble::column_to_rownames(var="Feature.ID") %>%
  tidyr::separate(Taxon, into=c( "Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
tax <- tax[row.names(dat),] %>%
  as.matrix() %>%
  tax_table()

tree <- phyloseq::read_tree("results/phylogenetic_tree/tree.nwk")

samp <- openxlsx::read.xlsx("inst/extdata/DN19196_191024_M02212_0230_000000000-CN2NW.xlsx",
                           sheet=1, startRow=2, rowNames=TRUE)[c( "Sample.Name", "Sample.Treatment")] 
samp$Number <-   sub(".*_", "", samp$Sample.Name)
names(samp) <- sub("Sample\\.", "", names(samp))
samp <- sample_data(samp[colnames(dat),])

physeq <- phyloseq::phyloseq(dat, tax, tree, samp)

usethis::use_data(physeq)
