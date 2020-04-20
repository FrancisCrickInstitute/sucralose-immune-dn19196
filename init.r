library(tidyverse)
library(dada2)
library(phyloseq)
library(usethis)
library(openxlsx)

conv_table <- openxlsx::read.xlsx("inst/extdata/conversion_table.xlsx", sheet=1)
all_samp <- samp <- openxlsx::read.xlsx("inst/extdata/DN19196_200125_M02212_0247_000000000-CWBLR.xlsx",
                                      sheet=1, startRow=2, rowNames=TRUE)[c( "Sample.Name", "Sample.Treatment", "Sample.Replicate.Group")] %>%
  tidyr::extract(Sample.Name, into=c("Week", "Batch"), ".*week(2|12).*_([0-9]{3})", remove=FALSE)
all_samp$Treatment <- with(all_samp, sub(" ", "_", tolower(ifelse(Sample.Treatment != "text", Sample.Treatment, Sample.Replicate.Group))))
all_samp <- all_samp[c("Week","Batch","Treatment","Sample.Name")]

conv_table <- cbind(conv_table, all_samp[match(conv_table$`2.week.sample.name`, all_samp$Sample.Name),])
ind <- dplyr::coalesce(
  match(all_samp$Sample.Name, conv_table$matching.week.12.name),
  match(all_samp$Sample.Name, conv_table$`2.week.sample.name`)
)
all_samp$Batch <- conv_table$Batch[ind]
all_samp$Mouse <- conv_table$mouse[ind]



physeqs <- list()

dsets <- c("results_1", "results")
dat <- list()
tax <- list()
tree <- list()

for (i in dsets) {
  
  dat[[i]] <- read.table(paste0(i, "/abundance_table/filtered/feature-table.tsv"), header=TRUE, comment.char="", skip=1, sep="\t") %>%
    tibble::column_to_rownames(var="X.OTU.ID") %>%
    as.matrix() %>%
    otu_table(taxa_are_rows=TRUE)
  
  
  tax[[i]] <- read.table(paste0(i, "/taxonomy/taxonomy.tsv"), header=TRUE, comment.char="", sep="\t") %>%
    dplyr::select(-Confidence) %>%
    dplyr::mutate(Taxon=gsub("D_[0-6]__","", Taxon)) %>%
    tibble::column_to_rownames(var="Feature.ID") %>%
    tidyr::separate(Taxon, into=c( "Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
  tax[[i]] <- tax[[i]][row.names(dat[[i]]),] %>%
    as.matrix() %>%
    tax_table()
  
  tree[[i]] <- phyloseq::read_tree(paste0(i, "/phylogenetic_tree/tree.nwk"))
  
  samp <- sample_data(all_samp[colnames(dat[[i]]),])
  
  physeqs[[i]] <- phyloseq::phyloseq(dat[[i]], tax[[i]], samp)
}


physeq <- merge_phyloseq(physeqs[[1]], physeqs[[2]])

usethis::use_data(physeq, overwrite=TRUE)
