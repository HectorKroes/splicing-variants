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

    if file -b --mime-type ${input_vcf} | grep -q x-gzip; then     
        zcat ${input_vcf} | grep -q "chr" && zcat ${input_vcf} | bcftools annotate --rename-chrs ${chr_format} -Oz -o ${input_vcf}
    elif file -b --mime-type ${input_vcf} | grep -q gzip; then
        zcat ${input_vcf} | grep -q "chr" && zcat ${input_vcf} | bcftools annotate --rename-chrs ${chr_format} -Oz -o ${input_vcf}
    else
        grep -q "chr" ${input_vcf} && bcftools annotate --rename-chrs ${chr_format} ${input_vcf} -Oz -o ${input_vcf}.gz
    fi

    bcftools index -f -t ${input_vcf.SimpleName}.vcf.gz
    """
}

process format_reference_files {

    /* Function to standardize reference files and provide
    a single copy of them to be used whenever necessary */ 
    
    stageInMode 'copy'
    stageOutMode 'copy'
    label 'spliceaiContainer'
    label 'inParallel'

    input:
        path fasta_file

    output:
        path("${fasta_file.SimpleName}.fa")

    script:
    """
    if file -b --mime-type ${fasta_file} | grep -q gzip; then
        zcat ${fasta_file} | sed 's/>chr/>/g' > ${fasta_file.baseName}
    else
        sed -i 's/>chr/>/g' ${fasta_file} > ${fasta_file}
    fi
    """

}

process annotate_variants {

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
    python ${python_file} ${input_vcf} ${params.spliceai_cutoff} ${params.squirls_cutoff} ${params.annotation_mode}
    """
}

workflow preanalysis {

  // Pre-analysis subworkflow

  take:
    input_vcfs
    chr_format
    fasta_file
 
  main: 
    input_files = format_input_files(input_vcfs, chr_format)
    fasta_ref = format_reference_files(fasta_file)

  emit: 
    input_files
    fasta_ref

}

workflow postanalysis {

    // Post-analysis subworkflow

  take:
    results
    annotation_script

  main:
    annotate_variants(results, annotation_script)
    
}