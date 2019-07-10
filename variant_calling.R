# The following R script demonstrates the usage of VariantTools R package 
# which provides a mechanisim for calling variants.

# align the reads
library(gmapR)
param <- GsnapParam(TP53Genome(), unique_only = TRUE, molecule = "DNA")
extdata.dir <- system.file("extdata", package="VariantToolsData")
first.fastq <- dir(extdata.dir, "first.fastq", full.names=TRUE)
last.fastq <- dir(extdata.dir, "last.fastq", full.names=TRUE)
output <- gsnap(first.fastq[1], last.fastq[1], param)
bam <- as(output, "BamFile")

# tally the variants
library(VariantTools)
data(repeats, package = "VariantToolsData")
genome(repeats) <- genome(TP53Genome())
param <- TallyVariantsParam(TP53Genome(), mask = repeats)
tallies <- tallyVariants(bam, param)

# filters
calling.filters <- VariantCallingFilters()
post.filters <- VariantPostFilters()

# calling variants
variants <- callVariants(tallies, calling.filters, post.filters)

# Alternative allele frequencies
variants$altFraction <-altDepth(variants) / totalDepth(variants)
library(ggplot2)
qplot(altFraction, geom = "density", data = as.data.frame(variants))

# How well our calls recapitulate the genotypes from 1000G project
data(geno, package = "VariantToolsData")

# Merge the expected frequencies of each alt with the variant calls.
naToZero <- function(x) ifelse(is.na(x), 0L, x)
addExpectedFreqs <- function(x) {
expected.freq <- geno$expected.freq[match(x, geno)]
x$expected.freq <- naToZero(expected.freq)
x
}
variants <- addExpectedFreqs(variants)

# dbSNP concordance
vcfPath <- system.file("extdata", "dbsnp-p53.vcf.gz", package = "VariantToolsData")
param <- ScanVcfParam(fixed = "ALT", info = NA, geno = NA)
dbSNP <- as(readVcf(vcfPath, param, genome = "TP53_demo"), "VRanges")
dbSNP <- dbSNP[!isIndel(dbSNP)]
variants$dbSNP <- variants %in% dbSNP
xtabs(~ dbSNP + expected.freq, mcols(variants))