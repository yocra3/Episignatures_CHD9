#! Rscript
#'#################################################################################
#'#################################################################################
#' Run methylation QC data
#' This code performs quality control to methylation data. 
#' Input: 
#' - path to config file with the parameters for the QC
#' - Name of the dataset to name the output files
#' Important decisions:
#' - Remove samples based on QC
#' - Use values from meffil vignette in all parameters - see config file
#'#################################################################################
#'#################################################################################

## Load libraries ####
library(meffil)
library(minfi)
library(ggplot2)
library(dplyr)

## Set number of cores
options(mc.cores = 10)

outPrefix <- "results/preprocess/script/chd9"
## Prepare sample sheet ####
samplesheet <- meffil.create.samplesheet(path = "data/1.Raw.data/idat/")

### Load and adapt samples data
load("results/phenotypes/script_phenotypes.Rdata")

## Merge both (Check in other datasets)
comb$Sample_Name <- paste(comb$Sentrix.Barcode, comb$Sample.Section, sep = "_")
combSheet <- left_join(dplyr::select(samplesheet, -Sex), comb, by = "Sample_Name")
combSheet$Sex <- combSheet$Sexo

## Generate QC report
## Load methylation data and QC ####
qc.objects <- meffil.qc(combSheet, verbose = TRUE,   cell.type.reference = "blood gse35069 complete")
qc.parameters <- meffil.qc.parameters(
  beadnum.samples.threshold             = 0.1,
  detectionp.samples.threshold          = 0.1,
  detectionp.cpgs.threshold             = 0.1, 
  beadnum.cpgs.threshold                = 0.1,
  sex.outlier.sd                        = 5,
  snp.concordance.threshold             = 0.95,
  sample.genotype.concordance.threshold = 0.8
)

qc.summary <- meffil.qc.summary(
  qc.objects,
  parameters = qc.parameters
)

outs <- c()

## Define function to select outliers
filterOutliers <- function(outlier){
  subset(outlier, issue %in% c("Control probe (dye.bias)",
                               "Methylated vs Unmethylated",
                               "X-Y ratio outlier",
                               "Low bead numbers",
                               "Detection p-value",
                               "Sex mismatch",
                               "Genotype mismatch",
                               "Control probe (bisulfite1)",
                               "Control probe (bisulfite2)")
  )
}

## Remove bad samples based on QC report and rerun QC
outlier <- qc.summary$bad.samples
outlier <- filterOutliers(outlier)
round <- 1
save(qc.objects, file = paste0(outPrefix, ".qc.objects.round", round, ".Rdata"))
save(qc.summary, file = paste0(outPrefix, ".qcsummary.round", round, ".Rdata"))
meffil.qc.report(qc.summary, output.file = paste0(outPrefix, ".methylationQC.raw.html"))


while (nrow(outlier)> 0){
  outs <- rbind(outs, outlier)
  round <- round + 1
  qc.objects <- meffil.remove.samples(qc.objects, outlier$sample.name)
  save(qc.objects, file = paste0(outPrefix,".qc.objects.round", round, ".Rdata"))
  
  qc.summary <- meffil.qc.summary(qc.objects, parameters = qc.parameters)
  save(qc.summary, file = paste0(outPrefix, ".qcsummary.round", round, ".Rdata"))
  outlier <- qc.summary$bad.samples
  outlier <- filterOutliers(outlier)
}
save(qc.objects, file = paste0(outPrefix, ".qc.objects.clean.Rdata"))
save(qc.summary, file = paste0(outPrefix, ".qcsummary.clean.Rdata"))
meffil.qc.report(qc.summary, output.file = paste0(outPrefix, ".methylationQC.clean.html"))

## Report filtered samples and probes
write.table(outs, file = paste0(outPrefix, ".removed.samples.txt"), quote = FALSE, row.names = FALSE,
            sep = "\t")
write.table(qc.summary$bad.cpgs, file = paste0(outPrefix, ".removed.probes.txt"), quote = FALSE, row.names = FALSE,
            sep = "\t")

## Run functional normalization ####
## To be changed in other projects
### Select number PCs (run, see plot and adapt pcs number)
y <- meffil.plot.pc.fit(qc.objects)
ggsave(y$plot, filename = paste0(outPrefix, ".pc.fit.pdf"), height = 6, width = 6)

norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs = 3)

## Add predicted sex as sample sheet variable
for (i in seq_len(length(norm.objects))){
  norm.objects[[i]]$samplesheet$pred.sex <- norm.objects[[i]]$predicted.sex
}

save(norm.objects, file = paste0(outPrefix, ".norm.obj.pc.Rdata"))

norm.beta <- meffil.normalize.samples(norm.objects, cpglist.remove = qc.summary$bad.cpgs$name, verbose = TRUE)

## Check covariables
beta.pcs <- meffil.methylation.pcs(norm.beta, probe.range = 40000)

batch_var <- c("Slide", "sentrix_row",  "sentrix_col", "Grupo", "Grupo2", "Sex", 
               "Estatus", "Familia")


## Define normalization parameters
norm.parameters <- meffil.normalization.parameters(
  norm.objects,
  variables = batch_var,
  control.pcs = seq_len(8),
  batch.pcs = seq_len(8),
  batch.threshold = 0.01
)
norm.summary <- meffil.normalization.summary(norm.objects, pcs = beta.pcs, parameters = norm.parameters)
save(norm.summary, file = paste0(outPrefix, ".norm.summary.Rdata"))
meffil.normalization.report(norm.summary, output.file = paste0(outPrefix, ".methylationQC.normalization.html"))

## Create GenomicRatioSet
rownames(combSheet) <- combSheet$Sample_Name
combSheet.final <- combSheet[colnames(norm.beta), ]
## Add predicted sex
combSheet.final$pred.sex <- vapply(norm.objects, function(x) x$predicted.sex, character(1))
gset <- makeGenomicRatioSetFromMatrix(norm.beta, pData = combSheet.final,
                                      array = "IlluminaHumanMethylationEPIC",
                                      annotation = "ilm10b4.hg19")


dp <- meffil.load.detection.pvalues(qc.objects)
detP <- dp[rownames(gset), colnames(gset)]

gset <- gset[rowMeans(detP < 2e-16) > 0.95, ]
save(gset, file = paste0(outPrefix, ".normalizedRaw.GenomicRatioSet.Rdata"))

