#! /usr/local/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Get the CpGs more associated with each of the groups
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
gsetfile <- args[1]

## Load libraries ####
library(minfi)
library(limma)
library(tidyverse)

## Define functions
getFeatures <- function(gset, group){
  
  gset.filt <- gset[, gset$Grupo2 %in% c(group, "Control")]
  model <- model.matrix(~ Grupo2 + factor(Sentrix.Barcode), colData(gset.filt))
  
  lmfit <- lmFit(getBeta(gset.filt), design = model)
  lmFite <- eBayes(lmfit)
  tab <- topTable(lmFite, n = Inf, coef = 2)
  tab.fil <- subset(tab, abs(logFC) > 0.20 & P.Value < 1e-3)
  if (nrow(tab.fil) > 1000){
    tab.fil <- tab.fil[1:1000, ]
  }
  rownames(tab.fil)
}

# Load data
load(gsetfile)

## Run models
groups <- c("CHD9", "REPS2", "Caso")
names(groups) <- groups

features <- lapply(groups[groups %in% unique(gset$Grupo2)], getFeatures, gset = gset)
save(features, file = "group_features.Rdata")