#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Slivar subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Implements the usage of the Slivar v0.3.0 software to annotate 
    allele frequencies
----------------------------------------------------------------------------------------
*/

process annotate_allele_frequency {

    /* Uses Slivar to annotate Gnomad allele
    frequencies */
    
    label 'inParallel'

    input:
        each input_vcf
        path gnomad_annotation

    output:
        path "${input_vcf.baseName}.vcf"

    script:
    """
    ./slivar expr --js slivar-functions.js -g ${gnomad_annotation} --vcf ${input_vcf.baseName}.vcf > ${input_vcf.baseName}.vcf
    """
}

workflow slivar {

  // Slivar subworkflow

  take:
    input_vcf
    gnomad_annotation
 
  main: 
    slivar_annotation = annotate_allele_frequency(input_vcf, gnomad_annotation)

  emit: 
    slivar_annotation

}