#' ---
#' title: "16S Analysis for Fabio Zani"
#' author: "Gavin Kelly"
#' output:
#'   html_document2:
#'     toc: true
#'     code_folding: hide
#'     css: styles.css
#' ---
#'

#+ init,results='hide', warning=FALSE, error=FALSE, message=FALSE
devtools::load_all()
library(tidyverse)
library(ggrepel)
library(xlsx)
library(phyloseq)
library(DESeq2)


param <- ParamList$new(title="Fabio's 16s Analysis",
                      script="analyse.r")
set.seed(param$get("seed"))

caption <- captioner()
knitr::opts_chunk$set(warning=FALSE, error=FALSE, message=FALSE,
                      dev=c("png","pdf"), out.width="80%",
                      results='asis')

#' # Preparation
#'
#' We take the output from the [ampliseq nf-core pipeline](https://github.com/nf-core/ampliseq)
#' generate aggregated coumts per treatment and per batch.  We also create a copy of the original
#' data scaled to total counts per sample. We discard OTUs that are aren't universally very low.
#' 
#+ aggregation, fig.cap=caption()
data(physeq, package="babs16s")
physeq <- list(base=physeq)
physeq$treatment <- merge_samples(physeq$base, sample_data(physeq$base)$Treatment)
descrip$treatment <- "Aggregated to treatments"
physeq$batch <- merge_samples(physeq$base, sample_data(physeq$base)$Batch)
descrip$batch <- "Aggregated to batches"

physeq$base_prop <- transform_sample_counts(physeq$base, function(OTU) OTU/sum(OTU)) # %>%
descrip$treatment <- "Scaled within samples"
param$set("low_count", value=2, description="Low counts are &le;  {}")
param$set("min_n_sample", value=3, description="Count must be be non-low in at least {} samples")

physeq$base <-  filter_taxa(physeq$base, filterfun(kOverA(k=param$get("min_n_sample"),A=param$get("low_count"))), TRUE)

#' # Richness of samples per group
#'
#' There are a number of measures of &alpha;-diversity, here we look at a few measures, examining
#' any patterns they might exhibit between batches and treatments:
#' 
#+  alpha-diversity, fig.cap=caption()

plot_richness(physeq$base, measures=c("Shannon", "Simpson"),  x="Treatment", color="Batch")
caption("Grouped by Treatment")
plot_richness(physeq$base, measures=c("Shannon", "Simpson"),  color="Treatment", x="Batch")
caption("Grouped by Batch")

#' # Taxon proportions {.tabset}
#'
#' We look at the majority OTUs within each taxon in turn, looking to see if there are any
#' big changes between samples, treatment groups and batches.
#'
#+ taxon-prop, fig.cap=caption()
param$set("topk", 5, "Each sample contributes {} OTUs")


stack_by <- function(agg, by, topk=5) {
  isblank <- is.na(as.matrix(tax_table(agg))[,by])
  tax_table(agg)[isblank,by] <- ""
  agg <- merge_taxa(agg, which(isblank), 1)
  agg <- tax_glom(agg, by, NArm=FALSE)
  istop <- genefilter_sample(agg, filterfun_sample(topk(topk)), 1)
  agg <- merge_taxa(agg, which(!istop), 1)
  tax_table(agg)[which(!istop)[1], by] <- "Other"
  print(plot_bar(agg, fill=by))
}


for (i in rank_names(physeq$base)[-1]) {
  cat("\n\n## ", i, "\n", sep="")
  stack_by(physeq$base, by=i, topk=param$get("topk"))
  caption(paste0("Major ", i, " in raw samples"))
  stack_by(physeq$base_prop, by=i, topk=param$get("topk"))
  caption(paste0("Major ", i, " in scaled samples"))
  stack_by(physeq$treatment, by=i, topk=param$get("topk"))
  caption(paste0("Major ", i, " in treatment groups"))
  stack_by(physeq$batch, by=i, topk=param$get("topk")) 
  caption(paste0("Major ", i, " in batches"))
}

#' # Ordination
#'
#' We want to be able to visualise the separation between samples. There are many
#' different meanings of 'distance between two samples' and here we examine a few of
#' them, and try to represent them as faithfully as possible in a two-dimensional figure.
#' 
#+ ordination, fig.cap=caption()

ord_nmds_bray <- ordinate(physeq$base_prop, method="NMDS", distance="bray", verbose=FALSE, trace=FALSE)
pl <- plot_ordination(physeq$base_prop, ord_nmds_bray,  title="Bray NMDS", color="Treatment") +
  coord_fixed()
gg <- ggplot_build(pl)$layout$panel_params
nmajor <- max(length(gg[[1]]$x.major_source), length(gg[[1]]$y.major.source))

pl <- pl + scale_x_continuous(breaks=scales::pretty_breaks(n=nmajor), labels=NULL) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=nmajor), labels=NULL) +
  labs(x="", y="") + theme_light() +
  theme(axis.ticks=element_blank(),
        plot.title = element_text(hjust = 0.5)
        )

print(pl + aes(colour=Batch, shape=Treatment))
caption("Bray distance, with an NMDS representation")


#' # Differential analysis
#'
#' As an exmaple, we can start to find differential OTUs. We assume
#' there is an underlying batch effect, and look to see where there is
#' a treatment effect once that batch effect has been accounted for.
#' We can test for any specific pairwise differences between
#' treatments, but for now we're going to look to see which OTUs
#' exhibit some significant variability across the four treatments
#' as a whole (ie the null hypothesis being they're all the same).
#'
#' Here we plot the fold-change between water and glucose
#' (arbitrary) for the statistically significant OTUs, showing which
#' phyla they belong to, labelling them by their genus.
#' 
#+ deseq, fig.cap=caption()
param$set("alpha", 0.05)
dds <- phyloseq_to_deseq2(physeq$base, ~  Batch + Treatment)
tt <- tax_table(physeq$base)
mcols(dds) <- tt
dds <- list(base=dds)
descrip <- list(base="Unscaled samples")


dds$base <- DESeq(dds$base, sfType="ratio", test="LRT", reduced = ~ Batch)
mcols(dds$base) <- cbind(mcols(dds$base), results(dds$base, alpha=param$get("alpha")))


diff_genus <- function(obj, alpha=0.05) {
  as.data.frame(mcols(obj)) %>%
    dplyr::filter(padj<alpha) %>%
    mutate(Phylum = fct_reorder(Phylum, log2FoldChange,max),
           Genus = fct_reorder(Genus, log2FoldChange, max)) %>%
    ggplot(aes(x=Phylum, y=log2FoldChange, label=Genus)) +
    geom_point(size=2) +
    geom_text_repel() +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
}

genus_scale <- setNames(RColorBrewer::brewer.pal(dplyr::n_distinct(mcols(dds$base)$Genus), "Set1"),
                        names(sort(table(mcols(dds$base)$Genus), decreasing=TRUE)))

stats <- map(dds, ~ table(sign(mcols(.)$log2FoldChange), mcols(.)$padj<param$get("alpha")))

imap(dds, ~diff_genus(.x, alpha=param$get("alpha")) +
                     labs(title=descrip[.y],
                          subtitle=sprintf("%d up, %d down", stats[[.y]]["1", "TRUE"], stats[[.y]]["-1", "TRUE"]))
)
caption("Differential OTUs")

