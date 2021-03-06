## Analysis overview (methods)

The fastq files were processed using DADA2 (v 1.18, ref 1), truncating
the forward (respectively reverse) reads to 280 and 210 bases and
trimming them by 17 and 21 bases respectively, with a maximum of two
expected errors. Taxa and species assignment was carried out using
v132 of the silva database.

The processed data were then analysed in R 4.0.3 (ref 2) with the
phyloseq package (ref 3), aggregating the count data to the genus
level. DESeq2 (v1.30, ref 4) was used to estimate the log-fold changes
and p-values between experimental groups whilst accounting for an
observed batch effect that crossed the experimental groups.

### References

1. Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP
(2016). “DADA2: High-resolution sample inference from Illumina amplicon
data.” _Nature Methods_, *13*, 581-583. doi: 10.1038/nmeth.3869 (URL:
https://doi.org/10.1038/nmeth.3869).

2. R Core Team (2020). R: A language and environment for statistical
computing. R Foundation for Statistical Computing, Vienna, Austria.
URL https://www.R-project.org/.

3. phyloseq: An R package for reproducible interactive analysis and
graphics of microbiome census data. Paul J. McMurdie and Susan Holmes
(2013) PLoS ONE 8(4):e61217.

4. Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change
and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550
(2014)


## File description

`run_dada.r` does the preprocessing, expecting the fastq files (from
SRA) to be in the `fastq` subdirectory.

`init.r` brings the results of the preprocessing into an R object for
downstream analysis.

`analyse.r` does the per-mouse and per-experimental-condition analysis
(both descriptive, and also exploratory), summarises the data into
fold-changes and p-values for various contrasts.

`data.r` produces various figures for the publication.
