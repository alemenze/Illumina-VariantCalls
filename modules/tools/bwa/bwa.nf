#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process bwa_index {
    tag "${meta}"
    label 'process_medium'

    publishDir "${params.outdir}/bwa/",
        mode: "copy",
        overwrite: true,
        saveAs: { filename -> filename }

    container "alemenze/bwa-tools"

    input:
        path(genome)
    
    output:
        path('*'), emit: index

    script:
        """
        bwa index $genome
        mkdir Index && mv ${genome}* Index
        """

}

process bwa_align {
    tag "${meta}"
    label 'process_high'

    publishDir "${params.outdir}/bwa/",
        mode: "copy",
        overwrite: true,
        saveAs: { filename -> filename }

    container "alemenze/bwa-tools"

    input:
        tuple val(meta), path(reads)
        path(genome)
        path(index)
    
    output:
        tuple val(meta), path("*.sam"),     emit: aligned_sam

    script:
        """
        bwa mem -t ${task.cpus} $index/${genome} $reads > ${meta}.sam
        """

}

process sam_sort {
    tag "${meta}"
    label 'process_medium'

    publishDir "${params.outdir}/bwa/",
        mode: "copy",
        overwrite: true,
        saveAs: { filename -> filename }

    container "alemenze/bwa-tools"

    input:
        tuple val(meta), path(sam)
    
    output:
        tuple val(meta), path('*.sorted.{bam,bam.bai}'), emit: aligned_bams
        tuple val(meta), path('*.sorted.bam'), emit: aligned_bam
        tuple val(meta), path("*{flagstat,idxstats,stats}"),   emit: logs

    script:
        """
        samtools view -hSbo ${meta}.bam $sam
        samtools sort ${meta}.bam -o ${meta}.sorted.bam
        samtools index ${meta}.sorted.bam
        samtools flagstat ${meta}.sorted.bam > ${meta}.sorted.bam.flagstat
        samtools idxstats ${meta}.sorted.bam > ${meta}.sorted.bam.idxstats
        samtools stats ${meta}.sorted.bam > ${meta}.sorted.bam.stats
        """

}

process mpileup {
    tag "${meta}"
    label 'process_medium'

    publishDir "${params.outdir}/vcfs/",
        mode: "copy",
        overwrite: true,
        saveAs: { filename -> filename }

    container "alemenze/bwa-tools"

    input:
        tuple val(meta), path(bam) 
        path(index)
        path(genome)
        path(gtf)
    
    output:
        tuple val(meta), path('*.variants.vcf'), emit: vcf

    script:
        """
        bcftools mpileup -f $index/${genome} $bam | bcftools call -mv -Ov > variants_temp.vcf
        bedtools intersect -a $gtf -b variants_temp.vcf -wa -u > ${meta}.variants.vcf
        """
}