#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process bwa_index {
    tag "${fasta}"
    label 'process_medium'

    publishDir "${params.outdir}/bwaindex/",
        mode: "copy",
        overwrite: true,
        saveAs: { filename -> filename }

    container "alemenze/bwa-tools"

    input:
        path(genome)
    
    output:
        path('bwa'), emit: index

    script:
        """
        mkdir bwa
        bwa index $genome -p bwa/${genome.baseName}
        """

}

process bwa_align {
    tag "${meta}"
    label 'process_high'

    publishDir "${params.outdir}/bwa/${meta}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename -> filename }

    container "alemenze/bwa-tools"

    input:
        tuple val(meta), path(reads)
        path index
        path(genome)
        path(gtf)
    
    output:
        tuple val(meta), path("*.sam"),     emit: aligned_sam
        tuple val(meta), path('*.sorted.{bam,bam.bai}'), emit: aligned_bams
        tuple val(meta), path('*.sorted.bam'), emit: aligned_bam
        tuple val(meta), path("*{flagstat,idxstats,stats}"),   emit: logs
        tuple val(meta), path('*.variants.vcf'), emit: vcf

    script:
        """
        INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`
        bwa mem -t ${task.cpus} \$INDEX $reads > ${meta}.sam

        samtools view -hSbo ${meta}.bam ${meta}.sam
        samtools sort ${meta}.bam -o ${meta}.sorted.bam
        samtools index ${meta}.sorted.bam
        samtools flagstat ${meta}.sorted.bam > ${meta}.sorted.bam.flagstat
        samtools idxstats ${meta}.sorted.bam > ${meta}.sorted.bam.idxstats
        samtools stats ${meta}.sorted.bam > ${meta}.sorted.bam.stats

        bcftools mpileup -f $genome ${meta}.sorted.bam | bcftools call -mv -Ov > variants_temp.vcf
        bedtools intersect -a $gtf -b variants_temp.vcf -wa -u > ${meta}.variants.vcf
        """

}
