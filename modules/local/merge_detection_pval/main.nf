// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MERGE_DETECTION_PVAL {

    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container 'yocra3/episignature_chd9:1.0'

    input:
    tuple path(detp1), path(detp2)

    output:
    path("combined.detectionPvalues.Rdata"), emit: detp

    script:
    """
    mergeDetectionPval.R $detp1 $detp2
    """
}
