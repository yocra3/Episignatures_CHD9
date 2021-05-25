#! /usr/local/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Merge GenomicRatioSet
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
detP1 <- args[1]
detP2 <- args[2]

## Load libraries
library(minfi)

## Prepare GenomicRatioSet ####
load(detP1)
det1 <- detP

load(detP2)
det2 <- detP

cpgs <- intersect(rownames(det1), rownames(det2))
samps <- setdiff(colnames(det2), colnames(det1))

detP <- cbind(det1[cpgs, ], det2[cpgs, samps])
save(detP, file = "combined.detectionPvalues.Rdata")
