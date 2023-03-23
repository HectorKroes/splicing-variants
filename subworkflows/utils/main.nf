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
    in a standardized format */ 
    
    stageInMode 'copy'

    input:
        path input_vcf

    output:
        path("${input_vcf.SimpleName}.vcf.gz")

    script:
    """
    if file -b --mime-type ${input_vcf} | grep -q x-gzip; then     
        bcftools index -f -t ${input_vcf}
    elif file -b --mime-type ${input_vcf} | grep -q gzip; then
        gunzip ${input_vcf}
        bcftools view ${input_vcf} -Oz -o ${input_vcf}.gz
    else     
        bcftools view ${input_vcf} -Oz -o ${input_vcf}.gz
    fi
    """
}