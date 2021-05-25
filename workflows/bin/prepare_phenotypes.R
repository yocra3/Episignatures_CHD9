#! /usr/local/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Merge phenotype files for QC
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
age <- args[1]
phenos <- args[2]

## Load libraries ####
library(tidyverse)

agedf <- read.delim(age)
phenodf <- read.delim(phenos)

comb <- phenodf %>%
  mutate(SampleID = Sample.ID) %>%
  select(-Sample.ID, - Index) %>%
  left_join(agedf) %>%
  mutate(Grupo = ifelse(Clínica == "X", "Caso", 
                        ifelse(Control == "X", "Control", "Portador")),
         Grupo2 = ifelse(Grupo == "Caso" & Mutación != "", Mutación, Grupo),
         Estatus = ifelse(Comentarios == "", "Probando", Comentarios))
save(comb, file = "phenotypes.Rdata")
write.table(comb, col.names = TRUE, quote = FALSE, file = "phenotypes.txt", 
            row.names = FALSE)