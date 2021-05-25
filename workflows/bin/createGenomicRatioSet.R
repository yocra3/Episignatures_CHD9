#! /usr/local/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Create GenomicRatioSets from matrix of betas
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
raw <- args[1]
phenos <- args[2]
outPrefix <- args[3]

## Load libraries
library(minfi)
library(tidyverse)

## Prepare GenomicRatioSet ####
raw_df <- read_table2(raw)
load(phenos)

beta_mat <- select(raw_df, ends_with("AVG_Beta")) %>%
  data.matrix()
rownames(beta_mat) <- raw_df$TargetID
colnames(beta_mat) <- gsub(".AVG_Beta", "", colnames(beta_mat))

rownames(comb) <- comb$SampleID

gset <- makeGenomicRatioSetFromMatrix(beta_mat, pData = comb[colnames(beta_mat), ],
                                      array = "IlluminaHumanMethylationEPIC",
                                      annotation = "ilm10b4.hg19")

## Make matrix of detection p-values ####
detP <- select(raw_df, ends_with(".Detection_Pval")) %>%
  data.matrix()
rownames(detP) <- raw_df$TargetID
colnames(detP) <- gsub(".Detection_Pval", "", colnames(detP))
detP <- detP[rownames(gset), colnames(gset)]

## Filter CpGs with call rate < 95%
gset <- gset[rowMeans(detP < 2e-16) > 0.95, ]
save(gset, file = paste0(outPrefix, ".normalizedBeta.GenomiRatioSet.Rdata"))
save(detP, file = paste0(outPrefix, ".detectionPvalues.Rdata"))