#' ---
#' title: "16S Analysis for Fabio Zani"
#' author: "Gavin Kelly"
#' output:
#'   bookdown::html_document2:
#'     toc: true
#'     code_folding: hide
#'     css: styles.css
#' ---
#'

#+ init,results='hide', warning=FALSE, error=FALSE, message=FALSE
devtools::load_all()
library(tidyverse)
library(ggrepel)
library(genefilter)
library(xlsx)
library(phyloseq)
library(DESeq2)
library(emmeans)

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
sample_data(physeq)$group <- paste0(sample_data(physeq)$Treatment, ":", sample_data(physeq)$Week)

## Change cross-family genus to be unique
tt <- tax_table(physeq)
ind <- tt[,"Genus"] %in% names(which(apply(table( tt[,"Family"], tt[,"Genus"])!=0, 2, sum)!=1))
tt[ind, "Genus"] <- paste(tt[ind, "Family"], tt[ind,"Genus"])
tax_table(physeq) <- tt


agg_mice <- merge_samples(physeq, "group")
sample_data(agg_mice)$Treatment <- sample_data(physeq)[match(sample_names(agg_mice), sample_data(physeq)$group),"Treatment"]
agg_mice <- tax_glom(agg_mice, taxrank="Genus")

## Aggregate at each level
aggs <- lapply(setNames(rank_names(physeq), rank_names(physeq))[2:6],
              function(rnk) {
                obj <- tax_glom(physeq, taxrank=rnk)
                taxa_names(obj) <- tax_table(obj)[,rnk]
                obj
              }
              )

contrs <- list(
  "Within Week" = expand.grid(Week=c(2,12), Treatment=c("glucose", "low_sucralose", "high_sucralose")),
  "Within Treatment" = expand.grid(Week=c(12), Treatment=c("water", "glucose", "low_sucralose", "high_sucralose"))
)

contrs$`Within Week` <- contrs$`Within Week` %>%
  mutate(numer1=paste("Treatment", Treatment, "vs_water", sep="_"),
         numer2= ifelse(Week==2, NA, paste0("Treatment",Treatment,".Week12"))
         )
contrs$`Within Treatment` <- contrs$`Within Treatment` %>%
  mutate(numer1=paste("Week", Week, "vs_2", sep="_"),
         numer2= ifelse(Treatment=="water", NA, paste0("Treatment",Treatment,".Week12"))
         )

ddss <- lapply(aggs, function(phy) {
  dds <- phyloseq_to_deseq2(phy, ~  Batch + Treatment*Week)
  dds$Week <- factor(dds$Week, levels=c(2,12))
  dds$Treatment <- factor(dds$Treatment, levels=c("water", "glucose", "low_sucralose", "high_sucralose"))
  dds <- DESeq(dds, sfType="ratio")
  lapply(contrs, function(x) {
    this_res <- list()
    for (i in 1:nrow(x)) {
      this_res[[i]]  <- results(dds, contrast=list(c(x$numer1[i], na.omit(x$numer2[i]))))
    }
    x$res <- this_res
    x
  }
  )
})

df <- psmelt(agg_mice) %>%
  select(-Sample.Name, -group, -Batch, -Mouse)

#write_csv(df, path="grouped_genus.csv")
#saveRDS(df, file="Grouped_genus.rds")

#ind <- df$Genus %in% names(which(apply(table( df$Family, df$Genus)!=0, 2, sum)!=1))
#df$Genus[ind] <- df$Family[ind]


## Put in a nested order of total abundance
df <- df %>%
  group_by(Phylum) %>% mutate(phylum=sum(Abundance)) %>%
  group_by(Class) %>% mutate(class=sum(Abundance)) %>%
  group_by(Order) %>% mutate(order=sum(Abundance)) %>%
  group_by(Family) %>% mutate(family=sum(Abundance)) %>%
  group_by(Genus) %>% mutate(genus=sum(Abundance)) %>%
  filter(row_number()==1) %>%
  ungroup() %>%
  arrange(desc(phylum),
          desc(class),
          desc(order),
          desc(family),
          desc(genus))

df$Phylum <- factor(df$Phylum, levels=unique(df$Phylum))
df$Class <- factor(df$Class, levels=unique(df$Class))
df$Order <- factor(df$Order, levels=unique(df$Order))
df$Family <- factor(df$Family, levels=unique(df$Family))
df$Genus <- factor(df$Genus, levels=unique(df$Genus))


#df <- df %>% pivot_wider(names_from=Week, values_from=Abundance, id_cols=-Sample)


hier <- list(
  Phylum = df %>%
    group_by(Phylum) %>%
    summarise(Abundance=sum(Abundance), rank="Phylum") %>%
    arrange(Phylum) %>% mutate(running_Abundance=cumsum(Abundance)) %>%
    mutate(level=as.character(Phylum)) %>%
    select(-Phylum),
  Class = df %>%
    group_by( Class) %>%
    summarise(Abundance=sum(Abundance), rank="Class") %>%
    arrange(Class) %>% mutate(running_Abundance=cumsum(Abundance)) %>%
    mutate(level=as.character(Class)) %>%
    select(-Class),
  Order = df %>%
    group_by( Order) %>%
    summarise(Abundance=sum(Abundance), rank="Order") %>%
    arrange(Order) %>% mutate(running_Abundance=cumsum(Abundance)) %>%
    mutate(level=as.character(Order)) %>%
    select(-Order),
  Family = df %>%
    group_by( Family) %>%
    summarise(Abundance=sum(Abundance), rank="Family") %>%
    arrange(Family) %>% mutate(running_Abundance=cumsum(Abundance)) %>%
    mutate(level=as.character(Family)) %>%
    select(-Family),
  Genus = df %>%
    group_by( Genus) %>%
    summarise(Abundance=sum(Abundance), rank="Genus") %>%
    arrange(Genus) %>% mutate(running_Abundance=cumsum(Abundance)) %>%
    mutate(level=as.character(Genus)) %>%
    select(-Genus)
  )

################################################################
#### Within-treatment (12 vs 2 hour)

add_lfc <- function(hierarchy, ddsList, ind) {
  my_hier <- list()
  for (lev in names(hierarchy)) {
  these <- list()
  dds <- ddsList[[lev]][[ind]]
  for (j in 1:nrow(ddsList[[lev]][[ind]])) {
    res <- dds$res[[j]][hierarchy[[lev]]$level,]
    these[[j]] <- cbind(hierarchy[[lev]],
                       lfc=res$log2FoldChange,
                       Treatment=as.character(dds$Treatment[j]),
                       Week=as.character(dds$Week[j])
                       )
  }
  my_hier[[lev]] <- do.call(rbind, these)
  }
  do.call(rbind, my_hier) %>%
    mutate(rank=factor(rank, levels=names(hierarchy)))
}

per_treatment_df <- add_lfc(hier, ddss, "Within Treatment")

ggplot(per_treatment_df,
       aes(xmin=running_Abundance-Abundance, xmax=running_Abundance,
           ymin=as.integer(rank), ymax=as.integer(rank)+0.75, fill=lfc)) +
  geom_rect() +
  facet_wrap(~Treatment) + 
  scale_fill_gradient2(name="LFC 12 vs 2",limits=c(-5,5), oob=scales::squish) +
  scale_y_continuous(name=NULL,breaks=seq(1.375, 5.375,1), labels=levels(per_treatment_df$rank)) + 
  theme_bw() + theme(panel.grid=element_blank())
  
  

################################################################
#### Per week ( vs water)

per_week_df <- add_lfc(hier, ddss, "Within Week")

ggplot(per_week_df,
       aes(xmin=running_Abundance-Abundance, xmax=running_Abundance,
           ymin=as.integer(rank), ymax=as.integer(rank)+0.75, fill=lfc)) +
  geom_rect() +
  facet_wrap(Week~Treatment) + 
  scale_fill_gradient2(name="LFC vs water",limits=c(-5,5), oob=scales::squish) +
  scale_y_continuous(name=NULL,breaks=seq(1.375, 5.375,1), labels=levels(per_week_df$rank)) + 
  theme_bw() + theme(panel.grid=element_blank())
  
################################################################


physeq <- list(base=physeq, species=tax_glom(physeq, taxrank="Species"))
descrip <- list(base="All samples separately", species="Compress equal species")


physeq$base_prop <- transform_sample_counts(physeq$base, function(OTU) OTU/sum(OTU)) # %>%
descrip$base_prop <- "Scaled within samples"
physeq$species_prop <- transform_sample_counts(physeq$species, function(OTU) OTU/sum(OTU)) # %>%
descrip$base_prop <- "Scaled species-level samples"


#' # Richness of samples per group
#'
#' There are a number of measures of &alpha;-diversity, here we look at a few measures, examining
#' any patterns they might exhibit between batches and treatments:
#' 
#+  alpha-diversity, fig.cap=caption()


plot_richness(physeq$base, measures=c("Shannon", "Simpson"),  x="Treatment", color="Week")
caption("Grouped by Treatment")
plot_richness(physeq$base, measures=c("Shannon", "Simpson"),  color="Treatment", x="Week")
caption("Grouped by Week")

#' # Taxon proportions {.tabset}
#'
#' We look at the majority OTUs within each taxon in turn, looking to see if there are any
#' big changes between samples, treatment groups and batches.
#'
#+ taxon-prop, fig.cap=caption()
param$set("low_count", value=2, description="Low counts are &le;  {}")
param$set("min_n_sample", value=3, description="Count must be be non-low in at least {} samples")

physeq$base <-  filter_taxa(physeq$base, filterfun(kOverA(k=param$get("min_n_sample"),A=param$get("low_count"))), TRUE)
physeq$species <-  filter_taxa(physeq$species, filterfun(kOverA(k=param$get("min_n_sample"),A=param$get("low_count"))), TRUE)
param$set("topk", 5, "Each sample contributes {} OTUs")


stack_by <- function(agg, by, topk=5, x="Mouse") {
  isblank <- is.na(as.matrix(tax_table(agg))[,by])
  tax_table(agg)[isblank,by] <- ""
  agg <- merge_taxa(agg, which(isblank), 1)
  agg <- tax_glom(agg, by, NArm=FALSE)
  istop <- genefilter_sample(agg, filterfun_sample(topk(topk)), 1)
  agg <- merge_taxa(agg, which(!istop), 1)
  tax_table(agg)[which(!istop)[1], by] <- "Other"
  invisible(plot_bar(agg, fill=by, x=s))
}




for (i in rank_names(physeq$base)[-1]) {
  cat("\n\n## ", i, "\n", sep="")
  print(stack_by(physeq$base_prop, by=i, topk=param$get("topk")) + facet_grid(Week~Treatment, scales="free_x"))
  caption(paste0("Major ", i, " in scaled samples"))
}

#' # Taxon Abundances {.tabset}
#'
#' Similarly, for absolute abundances:
#' 
#+ taxon-abun, fig.cap=caption()


for (i in rank_names(physeq$base)[-1]) {
  cat("\n\n## ", i, "\n", sep="")
  print(stack_by(physeq$base, by=i, topk=param$get("topk")) + facet_grid(Week~Treatment, scales="free_x"))
  caption(paste0("Major ", i, " in samples"))
}

#' # Select Genera {.tabset}
#'
#' For the given the list of genera, here are there abundances:
#'
#+ paper-genera, fig.cap=caption()

candidates <- list(
  list(Family="Ruminococcaceae", Genus="Ruminococcus 1"),
  list(Family="Lachnospiraceae", Genus="Roseburia"),
  list(Family="Staphylococcaceae", Genus="Staphylococcus"),
  list(Family="Clostridiaceae 1"),
  list(Family="Bacillaceae", Genus="Pseudogracilibacillus"),
  list(Family="Erysipelotrichaceae"),
  list(Family="Christensenellaceae")
)

for (i in candidates) {
  if ("Genus" %in% names(i)) {
    name <- paste( i$Family, i$Genus)
    ind <- subset_taxa(physeq$species, Family==i$Family & Genus==i$Genus)
    
  } else {
    name <- paste( i$Family)
    ind <- subset_taxa(physeq$species, Family==i$Family)
  }
  ind <- tax_glom(ind, "Genus")
  cat("\n\n## ", name, "\n", sep="")
  dat <- psmelt(ind) %>%
    dplyr::mutate(Mouse=as.factor(Mouse),
                  Week=factor(Week, levels=c("2","12"))
                  )

  pl <-   ggplot(dat, aes(x=Week, y=Abundance, group=Mouse, colour=Mouse)) +
    geom_point() + geom_line() +
    facet_grid(Genus~Treatment) +
    guides(colour=FALSE, group=FALSE)
  print(pl + stat_summary(aes(group=Treatment), fun.y=mean, geom="line", colour="black"))
  caption(name)

    pl <-   ggplot(dat, aes(x=Treatment, y=Abundance)) +
    geom_point() +
    facet_grid(Genus~Week) +
    guides(colour=FALSE, group=FALSE)
  print(pl + stat_summary(aes(group=Week), fun.y=mean, geom="line", colour="black"))
  caption(paste(name, " by week"))

}
             

#' # Ordination
#'
#' We want to be able to visualise the separation between samples. There are many
#' different meanings of 'distance between two samples' and here we examine a few of
#' them, and try to represent them as faithfully as possible in a two-dimensional figure.
#' 
#+ ordination, fig.cap=caption()

col_scheme <- c(water= "grey", 
               low_sucralose= "#ADD8E6",
               high_sucralose = "#00008B",
               glucose="#71bc78")

regress <- function(phy, model) {
  mat <- log(as.matrix(otu_table(phy))+0.1)
  df <- sample_data(phy)[colnames(mat),]
  for (i in 1:nrow(mat)) {
    df$y <- as.vector(mat[i,rownames(df)])
    fit <- lm(model, data=as.data.frame(unclass(df)))
    mat[i,] <- mat[i,]-t(predict(fit, type="term", terms="Batch"))
  }
  mat <- round(exp(mat)-0.1)
  otu_table(phy) <- otu_table(mat, taxa_are_rows = taxa_are_rows(phy))
  phy
}

to2 <- function(phy) {
  mat <- log(as.matrix(otu_table(phy))+0.1)
  df <- sample_data(phy)[colnames(mat),]
  is2 <- which(df$Week==2)
  ind2 <- is2[match(df$Mouse, df$Mouse[is2])]
  for (i in 1:nrow(mat)) {
    mat[i,] <- mat[i,]/mat[i,ind2]
  }
  mat <- round(exp(mat)-0.1)
  otu_table(phy) <- otu_table(mat, taxa_are_rows = taxa_are_rows(phy))
  subset_samples(phy, df$Week==12)
}
    
physeq$norm <- regress(physeq$species, model=y~Treatment+Week+Batch)
physeq$to2 <- to2(physeq$species)



ord <- list(method="MDS", distance="bray")
this_phy <- subset_samples(physeq$species)
this_phy <- subset_taxa(this_phy,  Genus!="Lachnospiraceae NK4A136 group")
ord_obj <- ordinate(this_phy, method=ord$method, distance=ord$distance)
pl <- plot_ordination(this_phy, ord_obj,  title=paste(ord, collapse=": "), color="Treatment") +
  coord_fixed()
gg <- ggplot_build(pl)$layout$panel_params
nmajor <- max(length(gg[[1]]$x.major_source), length(gg[[1]]$y.major.source))
pl <- pl + scale_x_continuous(breaks=scales::pretty_breaks(n=nmajor), labels=NULL) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=nmajor), labels=NULL) +
  labs(x="", y="") + theme_bw() + 
  theme(axis.ticks=element_blank(),
        plot.title = element_text(hjust = 0.5),
        text=element_text(family="TeXGyreHeros")
        ) + geom_point(size=3)
print(pl + facet_wrap(~Week) + scale_color_manual(values=col_scheme))
caption(paste0(ord$distance, "distance, with an ", ord$method, " representation"))





print(pl + aes(colour=Batch, shape=Treatment))
caption("Bray distance, with an NMDS representation, Batch pattern")



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

## genus_scale <- setNames(RColorBrewer::brewer.pal(dplyr::n_distinct(mcols(dds$base)$Genus), "Set1"),
##                         names(sort(table(mcols(dds$base)$Genus), decreasing=TRUE)))

stats <- map(dds, ~ table(sign(mcols(.)$log2FoldChange), mcols(.)$padj<param$get("alpha")))

imap(dds, ~diff_genus(.x, alpha=param$get("alpha")) +
                     labs(title=descrip[.y],
                          subtitle=sprintf("%d up, %d down", stats[[.y]]["1", "TRUE"], stats[[.y]]["-1", "TRUE"]))
)
caption("Differential OTUs")

