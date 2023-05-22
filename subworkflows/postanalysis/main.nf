#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Postanalysis subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Implements functions required in the post-analysis phase of the
    general pipeline, annotating the variants according their scores
----------------------------------------------------------------------------------------
*/

process annotate_variants {

    /* Uses python script to annotate final pipeline
    results into vcf files */
    
    publishDir params.o, mode: 'copy'
    label 'inParallel'

    input:
        each input_vcf
        path vcf_annotation

    output:
        path "splice_**.vcf"

    script:
    """
    python3 vcf_annotation.py ${input_vcf} ${params.spliceai_cutoff} ${params.squirls_cutoff} ${params.annotation_mode}
    """
}

workflow postanalysis {

    // Post-analysis subworkflow

  take:
    vcf_files
    python_file

  main:
    annotate_variants(vcf_files, python_file)
    
}