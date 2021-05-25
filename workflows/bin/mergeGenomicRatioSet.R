#! /usr/local/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Merge GenomicRatioSet
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
gr1 <- args[1]
gr2 <- args[2]

## Load libraries
library(minfi)

## Prepare GenomicRatioSet ####
load(gr1)
gset1 <- gset

load(gr2)
gset2 <- gset

cpgs <- intersect(rownames(gset1), rownames(gset2))
samps <- setdiff(colnames(gset2), colnames(gset1))

gset <- cbind(gset1[cpgs, ], gset2[cpgs, samps])
save(gset, file = "combined.normalizedBeta.GenomiRatioSet.Rdata")
