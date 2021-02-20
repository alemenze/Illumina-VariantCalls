#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition

process trimgalore {
    tag "${meta}"
    label 'process_medium'

    publishDir "${params.outdir}/trimgalore/",
        mode: "copy",
        overwrite: true,
        saveAs: { filename -> filename }

    container "quay.io/biocontainers/trim-galore:0.6.6--0"

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("*.fq.gz"),       emit: reads
        tuple val(meta), path("*report.txt"),   emit: log
        tuple val(meta), path("*.html"),        emit: html 
        tuple val(meta), path("*.zip") ,        emit: zip

    script:
        """
        trim_galore \\
            --cores ${task.cpus} \\
            --fastqc \\
            --paired \\
            --gzip \\
            $reads
        """
}