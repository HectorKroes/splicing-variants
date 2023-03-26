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
include { spliceai_annotate_precalculated_scores } from './subworkflows/spliceai'
include { spliceai_predict_de_novo_variants } from './subworkflows/spliceai'
include { spliceai_fuse_temporary_vcfs } from './subworkflows/spliceai'
include { squirls_predict_variant_effect } from './subworkflows/squirls'

// Defining input parameters

input_vcfs = Channel.
    fromPath(params.i)

indel_annotation = Channel
    .fromPath(params.indels, checkIfExists: true)

indel_annotation_index = Channel
    .fromPath(params.indels + '.tbi', checkIfExists: true)

snv_annotation = Channel
    .fromPath(params.snvs , checkIfExists: true)

snv_annotation_index = Channel
    .fromPath(params.snvs + '.tbi', checkIfExists: true)

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

    input = format_input_files(input_vcfs)

    pcs_channel = spliceai_annotate_precalculated_scores(format_input_files.out, indel_annotation, indel_annotation_index, snv_annotation, snv_annotation_index)

    spliceai_predict_de_novo_variants(pcs_channel.tbc_ch.combine(fasta_ref))

    temporary_vcfs = (pcs_channel.pcs_ch).join(spliceai_predict_de_novo_variants.out)

    spliceai_results = spliceai_fuse_temporary_vcfs(temporary_vcfs)

    squirls_predict_variant_effect(spliceai_results, squirls_db)

    filter_relevant_variants(squirls_predict_variant_effect.out, relevancy_filter_script)

}