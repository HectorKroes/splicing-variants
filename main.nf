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

include { preanalysis } from './subworkflows/utils'
include { spliceai } from './subworkflows/spliceai'
include { squirls } from './subworkflows/squirls'
include { postanalysis } from './subworkflows/utils'

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

fasta_file = Channel
    .fromPath(params.fa)

squirls_db = Channel
    .fromPath(params.sdb)

chr_format = Channel
    .fromPath('./internals/chr_dict.txt')

annotation_script = Channel
    .fromPath('./scripts/vcf_annotation.py')

results_folder = Channel
    .fromPath(params.o)

workflow {

    // Defining main workflow

    formatted_data = preanalysis(input_vcfs, chr_format, fasta_file)

    spliceai_results = spliceai(formatted_data.input_files, snv_annotation, indel_annotation, snv_annotation_index, indel_annotation_index, formatted_data.fasta_ref)

    squirls_results = squirls(spliceai.out, squirls_db)

    postanalysis(squirls_results, annotation_script)

}