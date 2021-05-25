// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FILTER_GENOMICRATIOSET {

    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignature_chd9:1.0'

    input:
    tuple val(prefix), path (gset), path (detP)
    path manifest

    output:
    path "*.Rdata", emit: results

    script:
    """
    filterGenomicRatioSet.R $gset $detP $manifest $prefix
    """
}
