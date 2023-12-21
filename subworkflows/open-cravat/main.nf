#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OpenCRAVAT subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Implements the usage of OpenCRAVAT v2.4.2 as a VCF annotator
----------------------------------------------------------------------------------------
*/

process open_cravat_annotations {
    
    /* Function receives input vcf, and necessary squirls files,
    annotates all the variants for which there are available
    scores and then saves them to the designated output folder*/

    label 'inSeries'

    input:
        path input_vcf

    output:
        path "output/*.vcf"

    script:
        """
        oc module install -y vest revel provean mutpred1 mutationtaster dbscsnv cadd cadd_exome aloft
        oc run ${input_vcf} -l hg38 -a vest revel provean mutpred1 mutationtaster dbscsnv cadd cadd_exome aloft  -d output/ -t vcf
        """

}

workflow open_cravat {

  // SQUIRLS subworkflow

  take:
    squirls_results
 
  main: 
    open_cravat_results = open_cravat_annotations(squirls_results)

  emit: 
    open_cravat_results

}