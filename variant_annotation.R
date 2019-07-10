# The following R script demonstrates the usage of VariantAnnotation R package  which is used for annotating variants.

# Install the required packages and libraries.
if (!requireNamespace("BiocManager", usquietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("VariantAnnotation")
BiocManager::install("SNPlocs.Hsapiens.dbSNP.20101109")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("PolyPhen.Hsapiens.dbSNP131")

# Load data of chromosome 22 from 1000 Genomes.
library(VariantAnnotation)
fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
vcf <- readVcf(fl)
rowRanges(vcf)
length(rowRanges(vcf))
geno(vcf)
header(vcf)
geno(vcf)$GT
rowRanges(vcf)[which( names(rowRanges(vcf))== "22:50307386_G/A"), ]
info(vcf)[1:5, 1:5]

# Compare quality measures between novel (i.e., not in dbSNP)
# and known (i.e., in dbSNP) variants and the variant type present in the file. 
library(SNPlocs.Hsapiens.dbSNP.20101109)
rd <- rowRanges(vcf)
seqlevels(rd) <- "ch22"
getSNPcount()
cRh22snps <- getSNPlocs("ch22")
dbch22snps <- sub("rs", "", names(rd)) %in% chf22snps$RefSNP_id
table(dbch22snps)
metrics <- data.frame(QUAL=qual(vcf), inDbSNP=dbch22snps,
                      AVT=info(vcf)$VT, LDAF=info(vcf)$LDAF, RSQ=info(vcf)$RSQ)
library(ggplot2)
ggplot(metrics, aes(x=RSQ, fill=inDbSNP)) +
  geom_density(alpha=0.5) +
  scale_x_continuous(name="MaCH / Thunder Imputation Quality") +
  scale_y_continuous(name="Density") +
  theme(legend.position="top")

# Select genomic coordinates
rng <- GRanges(seqnames="22", ranges=IRanges(
  start=c(50301422, 50
          989541),
  end=c(50312106, 51001328),
  names=c("gene_79087", "gene_644186")))
tab <- TabixFile(fl)
vcf_rng <- readVcf(tab, "hg19", param=rng)
head(rowRanges(vcf_rng), 9)
length(rowRanges(vcf))
length(rowRanges(vcf_rng))

# Locating variants (Genetic Regions)
library(BiocManager::TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(vcf) <- "chr22"
rd <- rowRanges(vcf)
length(rd)
loc <- locateVariants(rd, txdb, CodingVariants())
head(loc, 3)
length(loc)
allvar <- locateVariants(rd, txdb, AllVariants())
length(allvar)

# Did any coding variants match more than one gene?
splt <- split(mcols(loc)$GENEID, mcols(loc)$QUERYID)
table(sapply(splt, function(x) length(unique(x)) > 1))

# Summarize the number of coding variants by gene ID.
splt <- split(mcols(loc)$QUERYID, mcols(loc)$GENEID)
head(sapply(splt, function(x) length(unique(x))), 3)

# Find amino-acid coding changes
library(BSgenome.Hsapiens.UCSC.hg19)
coding <- predictCoding(vcf, txdb, seqSource=Hsapiens)
coding[5:7]
coding[mcols(coding)$CONSEQUENCE == "frameshift"]

# identify the non-synonymous variants and obtain the rsids.
nms <- names(coding)
idx <- mcols(coding)$CONSEQUENCE == "nonsynonymous"
nonsyn <- coding[idx]
names(nonsyn) <- nms[idx]
rsids <- unique(names(nonsyn)[grep("rs", names(nonsyn), fixed=TRUE)])

# Query the PolyPhen databaseccv
library(PolyPhen.Hsapiens.dbSNP131)
pp <- select(PolyPhen.Hsapiens.dbSNP131, keys=rsids,
             cols=c("TRAININGSET", "PREDICTION", "PPH2PROB"))
head(pp[!is.na(pp$PREDICTION), ])
unique(pp$PREDICTION)