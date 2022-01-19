#'#################################################################################
#'#################################################################################
# Nextflow commands
#'#################################################################################
#'#################################################################################
export TOWER_ACCESS_TOKEN=eyJ0aWQiOiAzMzE2fS5kNjdiNWM3OGVkYWUyYWE3ZjM1OWYxOGNlYTU2NDBhYmMzNjdlODY2

## Prepare genotype files for imputation
nextflow run workflows/runMethylationQC.nf --pheno_age data/Edadysexo.csv \
--inputTable data/input.tsv -profile docker -resume

## Define episignatures
nextflow run workflows/defineEpimutations.nf \
--gset results/preprocess/script/chd9.autosomic.filterAnnotatedProbes.GenomicRatioSet.Rdata \
-profile docker -resume
