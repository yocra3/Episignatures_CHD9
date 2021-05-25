// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CREATE_GENOMICRATIOSET {

    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignature_chd9:1.0'

    input:
    tuple val(prefix), path(raw), path (pheno)

    output:
    path("*.normalizedBeta.GenomiRatioSet.Rdata"), emit: gset
    path("*.detectionPvalues.Rdata"), emit: detp

    script:
    """
    createGenomicRatioSet.R $raw $pheno $prefix
    """
}
