#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process multiqc {
    tag "${meta}"
    label 'process_low'

    publishDir "${params.outdir}/multiqc",
        mode: "copy",
        overwrite: true,
        saveAs: { filename -> filename }

    container "quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"

    input:
        path ('fastqc/*')
        path ('trimgalore/fastqc/*')
        path ('trimgalore/*')
        path ('bwa/*')
        path ('kraken/*')
    
    output:
        path "*multiqc_report.html", emit: report
        path "*_data"              , emit: data

    script:
        """
        multiqc -f .
        """
}