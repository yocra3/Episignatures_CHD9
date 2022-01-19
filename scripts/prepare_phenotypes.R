#! /usr/local/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Merge phenotype files for QC
#'#################################################################################
#'#################################################################################

## Capture arguments
age <- "data/Edadysexo.csv"
pheno1 <- "data/1.Raw.data/Sample.Table.txt"
pheno2 <- "data/rawdata/Sample.Table.txt"

## Load libraries ####
library(tidyverse)

agedf <- read.delim(age)
phenodf1 <- read.delim(pheno1)[, -1]
phenodf2 <- read.delim(pheno2)[, -1]
phenodf <- rbind(phenodf1, phenodf2) %>%
  distinct() %>%
  filter(Sample.Group != "A20_0646")

## Manually correct phenotdata
phenodf$Sample.ID[phenodf$Sample.Group == "19_0215"] <- "19_0213" 
phenodf$Sample.ID[phenodf$Sample.Group == "19_0213"] <- "19_0215" 
phenodf <- subset(comb, Sample.ID != "11_0119") 

comb <- phenodf %>%
  mutate(SampleID = Sample.ID) %>%
  select(-Sample.ID) %>%
  left_join(agedf) %>%
  mutate(Grupo = ifelse(Clínica == "X", "Caso", 
                        ifelse(Control == "X", "Control", "Portador")),
         Grupo2 = ifelse(Grupo == "Caso" & Mutación != "", Mutación, Grupo),
         Estatus = ifelse(Comentarios == "", "Probando", Comentarios))
save(comb, file = "results/phenotypes/script_phenotypes.Rdata")
write.table(comb, col.names = TRUE, quote = FALSE, file = "results/phenotypes/script_phenotypes.txt", 
            row.names = FALSE)
