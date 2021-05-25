#!/usr/bin/env nextflow
/*
========================================================================================
Define episignatures
========================================================================================
*/

nextflow.enable.dsl = 2

////////////////////////////////////////////////////
/* --     Collect configuration parameters     -- */
////////////////////////////////////////////////////

ch_gset = file("${params.gset}")

date = java.time.LocalDate.now()

include { GET_GROUP_FEATURES } from '../modules/local/get_group_features/main.nf' addParams( options: [publish_dir: "episignatures/${date}"])

workflow  {

  GET_GROUP_FEATURES( ch_gset )

}
