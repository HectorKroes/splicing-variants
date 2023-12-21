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

include { preanalysis } from './subworkflows/preanalysis'
include { spliceai_cpu } from './subworkflows/spliceai-cpu'
include { spliceai_gpu } from './subworkflows/spliceai-gpu'
include { squirls } from './subworkflows/squirls'
include { open_cravat } from './subworkflows/open-cravat'
include { postanalysis } from './subworkflows/postanalysis'

// Defining input parameters

input_vcfs = Channel
    .fromPath(params.i, checkIfExists: true)

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

slivar_annotation = Channel
    .fromPath(params.af)

normalization_script = Channel
    .fromPath('./resources/usr/bin/vcf_normalization.py')

freq_filter_script = Channel
    .fromPath('./resources/usr/bin/vcf_freq_filter.py')

annotation_script = Channel
    .fromPath('./resources/usr/bin/vcf_annotation.py')

chr_format = Channel
    .fromPath('./resources/usr/src/chr_dict.txt')

results_folder = Channel
    .fromPath(params.o)

workflow {

    // Defining main workflow

    input_files = preanalysis(input_vcfs, chr_format, fasta_file, slivar_annotation, normalization_script, freq_filter_script)

    if ( params.gpu ) {

        spliceai_results = spliceai_gpu(input_files.formatted_vcfs, snv_annotation, indel_annotation, snv_annotation_index, indel_annotation_index, input_files.fasta_ref)

    } else {

        spliceai_results = spliceai_cpu(input_files.formatted_vcfs, snv_annotation, indel_annotation, snv_annotation_index, indel_annotation_index, input_files.fasta_ref)

    }

    squirls_results = squirls(spliceai_results, squirls_db)

    open_cravat_results = open_cravat(squirls_results)

    postanalysis(open_cravat_results, annotation_script)

}