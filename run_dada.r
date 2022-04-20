##' Initialize
library(dada2)
path <- "fastq"
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))

fnRs <- sub("_R1_", "_R2_", fnFs)

##' QC
sample.names <- sub(".*(ZAN.*)_S.*", "\\1", fnFs)
plotQualityProfile(fnFs[c(1:8)])
plotQualityProfile(c(fnFs[1], fnRs[1]))

##' Filter
path <- "scratch"
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,210),trimLeft=c(17,21),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)


##' Learn errors
errF <- learnErrors(filtFs, multithread=TRUE, verbose=TRUE, nbases=1e8)
errR <- learnErrors(filtRs, multithread=TRUE, verbose=TRUE)
plotErrors(errF, nominalQ=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

##' Sample Inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

## mergers <- vector("list", length(sample.names))
## names(mergers) <- sample.names
## for(sam in sample.names) {
##   cat("Processing:", sam, "\n")
##     derepF <- derepFastq(filtFs[[sam]])
##     ddF <- dada(derepF, err=errF, multithread=TRUE)
##     derepR <- derepFastq(filtRs[[sam]])
##     ddR <- dada(derepR, err=errR, multithread=TRUE)
##     merger <- mergePairs(ddF, derepF, ddR, derepR)
##     mergers[[sam]] <- merger
## }
## rm(derepF); rm(derepR)
##' Merge reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
save(merger, file="data/merger.rda")

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)


taxa <- assignTaxonomy(seqtab.nochim, "../data/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "../data/silva_species_assignment_v132.fa.gz")

library(DECIPHER)
dna <- DNAStringSet(getSequences(seqtab.nochim))
load("../data/SILVA_SSU_r132_March2018.RData")
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
        m <- match(ranks, x$rank)
        taxa <- x$taxon[m]
        taxa[startsWith(taxa, "unclassified_")] <- NA
        taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)
save(seqtab.nochim, taxa, taxid, file="../objects/dada_long.rdata")

