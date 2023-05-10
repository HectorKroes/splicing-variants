#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Main nextflow workflow for splicing variant analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This is a DSL 2 Nextflow workflow that implements SpliceAI and SQUIRLS 
    to identify genetic variants that potentially affect RNA splicing
----------------------------------------------------------------------------------------
*/

// Including processes from subworkflows

include { format_input_files } from './subworkflows/utils'
include { filter_relevant_variants } from './subworkflows/utils'
include { spliceai } from './subworkflows/spliceai'
include { squirls } from './subworkflows/squirls'

// Defining input parameters

input_vcfs = Channel
    .fromPath(params.i)

indel_annotation = Channel
    .fromPath(params.indels)

indel_annotation_index = Channel
    .fromPath(params.indels + '.tbi')

snv_annotation = Channel
    .fromPath(params.snvs)

snv_annotation_index = Channel
    .fromPath(params.snvs + '.tbi')

fasta_ref = Channel
    .fromPath(params.fa, checkIfExists: true)

squirls_db = Channel
    .fromPath(params.sdb, checkIfExists: true)

relevancy_filter_script = Channel
    .fromPath('./scripts/relevancy_filter.py')

results_folder = Channel
    .fromPath(params.o)

workflow {

    // Defining main workflow

    input_files = format_input_files(input_vcfs)

    spliceai_results = spliceai(input_files, snv_annotation, indel_annotation, snv_annotation_index, indel_annotation_index, fasta_ref)

    squirls_results = squirls(spliceai.out, squirls_db)

    annotation = filter_relevant_variants(squirls_results, relevancy_filter_script)

}