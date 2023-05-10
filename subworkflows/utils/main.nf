#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Utils subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Implements common utilities functions for the main workflow
----------------------------------------------------------------------------------------
*/

process format_input_files {

    /* Function to ensure input vcf file enters workflow
    in a standardized format: in bgzip compression, and
    with chr formatting compatible with Illumina files */ 
    
    stageInMode 'copy'
    label 'spliceaiContainer'
    label 'inParallel'

    input:
        path input_vcf
        path chr_format

    output:
        tuple path("${input_vcf.SimpleName}.vcf.gz"), path("${input_vcf.SimpleName}.vcf.gz.tbi")

    script:
    """

    if file -b --mime-type ${input_vcf} | grep -q "gzip"; then
        zcat ${input_vcf} | grep -q ">chr" && zcat ${input_vcf} | bcftools annotate --rename-chrs ${chr_format} -o ${input_vcf}
    else
        grep -q ">chr" ${input_vcf} && bcftools annotate --rename-chrs ${chr_format} ${input_vcf} -o ${input_vcf}
    fi

    if file -b --mime-type ${input_vcf} | grep -q x-gzip; then     
        true
    elif file -b --mime-type ${input_vcf} | grep -q gzip; then
        gunzip ${input_vcf}
        bcftools view ${input_vcf} -Oz -o ${input_vcf}.gz
    else     
        bcftools view ${input_vcf} -Oz -o ${input_vcf}.gz
    fi

    bcftools index -f -t ${input_vcf.SimpleName}.vcf.gz
    """
}

process filter_relevant_variants {

    /* Uses python script to annotate final pipeline
    results into vcf files */
    
    publishDir params.o, mode: 'copy'
    label 'spliceaiContainer'
    label 'inParallel'

    input:
        each input_vcf
        path python_file

    output:
        path "relevant_**.vcf"

    script:
    """
    python relevancy_filter.py ${input_vcf} ${params.spliceai_cutoff} ${params.squirls_cutoff} ${params.annotation_mode}
    """
}