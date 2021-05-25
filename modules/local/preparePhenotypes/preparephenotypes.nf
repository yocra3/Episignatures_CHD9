// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PREPARE_PHENOTYPES {

    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignature_chd9:1.0'

    input:
    path pheno_age
    tuple val(prefix), path (pheno_disease)

    output:
    tuple val(prefix), path ("*.Rdata"), emit: rdata
    path "*.txt", emit: txt

    script:
    """
    prepare_phenotypes.R $pheno_age $pheno_disease
    """
}
