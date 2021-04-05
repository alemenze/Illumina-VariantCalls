#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

////////////////////////////////////////////////////
/* --              Parameter setup             -- */
////////////////////////////////////////////////////

if (params.samplesheet) {file(params.samplesheet, checkIfExists: true)} else { exit 1, 'Samplesheet file not specified!'}
if (params.genome) {file(params.genome, checkIfExists: true)} else { exit 1, 'Genome fasta file not specified!'}
if (params.gtf) {file(params.gtf, checkIfExists: true)} else { exit 1, 'GTF file not specified!'}

Channel
    .fromPath(params.samplesheet)
    .splitCsv(header:true)
    .map{ row -> tuple(row.sample_id, file(row.read1), file(row.read2))}
    .set { sample_reads }
Channel
    .fromPath(params.genome)
    .set {genome}
Channel
    .fromPath(params.gtf)
    .set {gtf}

////////////////////////////////////////////////////
/* --              IMPORT MODULES              -- */
////////////////////////////////////////////////////

include { fastqc } from '../tools/fastqc/fastqc'
include { trimgalore } from '../tools/trimgalore/trimgalore'
include { multiqc } from '../tools/multiqc/multiqc'
include { bwa_index } from '../tools/bwa/bwa'
include { bwa_align } from '../tools/bwa/bwa'
include { sam_sort } from '../tools/bwa/bwa'
include { mpileup } from '../tools/bwa/bwa'
include { Kraken } from '../subworkflows/kraken'

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow Alignment {
    
    fastqc(
        sample_reads
    )

    trimgalore(
        sample_reads
    )

    bwa_index(
        genome
    )

    bwa_align(
        trimgalore.out.reads,
        bwa_index.out.index.collect(),
        genome.collect(),
        gtf.collect()
    )
    
    Kraken(
        trimgalore.out.reads,
        'Kraken'
    )

    multiqc(
        fastqc.out.zip.collect{ it[1] },
        trimgalore.out.zip.collect{ it[1] },
        trimgalore.out.log.collect{ it[1] },
        bwa_align.out.logs.collect{ it[1] },
        Kraken.out.krakenreport.collect{ it[1] }
    )

} 
