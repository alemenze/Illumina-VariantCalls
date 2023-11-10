#!/usr/bin/env nextflow
/*
                              (╯°□°)╯︵ ┻━┻

========================================================================================
                            Variant calling workflow
========================================================================================
                  https://github.com/alemenze/Illumina-MTB-VariantCalls
*/

nextflow.enable.dsl = 2

def helpMessage(){
    log.info"""

    Usage:
    
    --samplesheet           Path to sample sheet
    --genome                Path to FASTA genome file
    --gtf                   Path to GTF annotation file
    
    """

}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

////////////////////////////////////////////////////
/* --              Parameter setup             -- */
////////////////////////////////////////////////////

////////////////////////////////////////////////////
/* --              IMPORT MODULES              -- */
////////////////////////////////////////////////////
include { Alignment } from './modules/main_workflows/alignment'
////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Full workflow demultiplexing through Trycycler
workflow {
    Alignment()
}
