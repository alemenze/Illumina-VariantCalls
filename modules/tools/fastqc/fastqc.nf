#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process fastqc {
    tag "${meta}"
    label 'process_low'

    publishDir "${params.outdir}/fastqc",
        mode: "copy",
        overwrite: true,
        saveAs: { filename -> filename }

    container "biocontainers/fastqc:v0.11.9_cv7"

    input:
        tuple val(meta), path(read1)
    
    output:
        tuple val(meta), path("*.zip"), emit: zip
        tuple val(meta), path("*.html"), emit: html

    script:
        """
        fastqc -q $read1
        """

}