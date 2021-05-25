#!/usr/bin/env nextflow
/*
========================================================================================
Run methylation QC of Episignature project with meffil
========================================================================================
*/

nextflow.enable.dsl = 2

////////////////////////////////////////////////////
/* --     Collect configuration parameters     -- */
////////////////////////////////////////////////////

ch_age = file("${params.pheno_age}")

manifest = file("${params.epic_manifest}")

Channel
        .fromPath(params.inputTable)
        .splitCsv(sep: "\t")
        .map { row -> [ row[0], file(row[1], checkIfExists: true)]}
        .set { ch_pheno }

Channel
  .fromPath(params.inputTable)
  .splitCsv(sep: "\t")
  .map { row -> [ row[0], file(row[2], checkIfExists: true)]}
  .set { ch_raw }

date = java.time.LocalDate.now()

include { PREPARE_PHENOTYPES } from '../modules/local/preparePhenotypes/preparephenotypes.nf' addParams( options: [publish_dir: "phenotypes/${date}"])
include { CREATE_GENOMICRATIOSET } from '../modules/local/creategenomicratioset/main.nf' addParams( options: [publish_dir: "preprocess/${date}", publish_files : ['GenomiRatioSet.Rdata':'']])
include { FILTER_GENOMICRATIOSET } from '../modules/local/filtergenomicratioset/main.nf' addParams( options: [publish_dir: "preprocess/${date}"])
include { MERGE_GENOMICRATIOSET } from '../modules/local/merge_genomicratioset/main.nf' addParams( options: [publish_dir: "preprocess/${date}"])
include { MERGE_DETECTION_PVAL } from '../modules/local/merge_detection_pval/main.nf' addParams( options: [publish_dir: "preprocess/${date}"])

workflow  {

  PREPARE_PHENOTYPES( ch_age, ch_pheno)

  ch_raw_input = ch_raw.join(PREPARE_PHENOTYPES.out.rdata)

  CREATE_GENOMICRATIOSET( ch_raw_input )
  MERGE_GENOMICRATIOSET( CREATE_GENOMICRATIOSET.out.gset.collect() )
  MERGE_DETECTION_PVAL( CREATE_GENOMICRATIOSET.out.detp.collect() )

  ch_grset = Channel.of('combined').concat(MERGE_GENOMICRATIOSET.out.gset).concat(MERGE_DETECTION_PVAL.out.detp)

  FILTER_GENOMICRATIOSET( ch_grset.collect(), manifest )

}
