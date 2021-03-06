#! /usr/local/bin/Rscript
#'#################################################################################
#'#################################################################################
#' This code filters GenomicRatioSet for analysis, by removing different types of
#' probes:
#' - Probes not measuring methylation (SNPs, CH probes)
#' - Crosshibridizing probes
#' - Probes with SNPs
#' - Probes in sexual chromosomes
#' 
#' For crosshibridizing probes and probes with SNPs, we used annotation from PMID: 27924034
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
gsetfile <-  "results/preprocess/script/chd9.normalizedRaw.GenomicRatioSet.Rdata"
qc_path <- "results/preprocess/script/chd9.qc.objects.clean.Rdata"
manifest <- "data/EPIC.hg19.manifest.rds"
outPrefix <- "results/preprocess/script/chd9"

## Load libraries ####
library(minfi)
library(meffil)

## Load dataset ####
load(gsetfile)
load(qc_path)
grAnnot <- readRDS(manifest)
ori <- gset

### Probes not measuring methylation
gset <- dropMethylationLoci(gset)
save(gset, file = paste0(outPrefix, ".allCpGs.GenomicRatioSet.Rdata"))

### Remove crosshibridizing and probes with SNPs
gset <- gset[!grAnnot[rownames(gset)]$MASK_general, ]
save(gset, file =  paste0(outPrefix, ".filterAnnotatedProbes.GenomicRatioSet.Rdata"))

### Remove probes in sexual chromosomes
gset <- gset[!seqnames(rowRanges(gset)) %in% c("chrX", "chrY"), ]
save(gset, file = paste0(outPrefix, ".autosomic.filterAnnotatedProbes.GenomicRatioSet.Rdata"))

## Create Initial and final dataset with missings
final <- gset
gset <- ori

detP <- meffil.load.detection.pvalues(qc.objects)
dp.f <- detP[rownames(gset), colnames(gset)]

beta <- getBeta(gset)
beta[dp.f > 2e-16] <- NA
assay(gset) <- beta
save(gset, file = paste0(outPrefix, ".allCpGs.withNA.GenomicRatioSet.Rdata"))

gset <- gset[rownames(final), colnames(final)]
save(gset, file = paste0(outPrefix, ".autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata"))

