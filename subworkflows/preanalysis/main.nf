#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Preanalysis subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Implements functions required in the pre-analysis phase of the
    general pipeline, mostly formatting-related
----------------------------------------------------------------------------------------
*/

process format_vcf_files {

    /* Function to ensure input vcf file enters workflow
    in a standardized format: in bgzip compression, and
    with chr formatting compatible with Illumina files */ 
    
    stageInMode 'copy'
    label 'inParallel'

    input:
        tuple path(input_vcf), path(chr_format)

    output:
        path("${input_vcf.SimpleName}.vcf.gz")

    script:
    """
    if file -b --mime-type ${input_vcf} | grep -q x-gzip; then     
        zcat ${input_vcf} | grep -q "chr" && zcat ${input_vcf} | bcftools annotate --rename-chrs ${chr_format} -Oz -o ${input_vcf}
    elif file -b --mime-type ${input_vcf} | grep -q gzip; then
        zcat ${input_vcf} | grep -q "chr" && zcat ${input_vcf} | bcftools annotate --rename-chrs ${chr_format} -Oz -o ${input_vcf}
    elif grep -q "chr" ${input_vcf}; then
        bcftools annotate --rename-chrs ${chr_format} ${input_vcf} -Oz -o ${input_vcf}.gz
    else
        bcftools view ${input_vcf} -Oz -o ${input_vcf}.gz
    fi
    """
}

process format_reference_files {

    /* Function to standardize reference files and provide
    a single copy of them to be used whenever necessary */ 
    
    label 'inParallel'

    input:
        path fasta_file

    output:
        path "${fasta_file.SimpleName}.fa", emit: fa
        path "${fasta_file.SimpleName}.fa.fai", emit: fai

    script:
    """
    zcat ${fasta_file} | sed 's/>chr/>/g' > ${fasta_file.SimpleName}.fa
    samtools faidx ${fasta_file.SimpleName}.fa
    """

}

process normalize_files {

    /* Function to normalize and index vcfs using bcftools */ 

    input:
        each vcf_file
        path fasta_file
        path normalization_script

    output:
        path "${vcf_file.SimpleName}.vcf.gz"

    script:
    """
    bcftools norm -cs -f ${fasta_file} ${vcf_file} -Oz -o normalized_vcf.vcf.gz
    gunzip normalized_vcf.vcf.gz 
    python3 ${normalization_script} normalized_vcf.vcf
    bcftools view normalized_vcf.vcf -Oz -o ${vcf_file.SimpleName}.vcf.gz
    """

}

process slivar_annotation {

    /* Uses Slivar to annotate Gnomad allele
    frequencies */
    
    stageInMode 'copy'
    label 'inParallel'

    input:
        tuple path(input_vcf), path(gnomad_annotation), path(fasta_bgzip), path(freq_filter_script)

    output:
        path "${input_vcf.SimpleName}.vcf.gz"
    script:
    """
    gunzip ${input_vcf}
    slivar expr -g ${gnomad_annotation} --vcf ${input_vcf.SimpleName}.vcf > slivar_vcf.vcf
    python3 vcf_freq_filter.py slivar_vcf.vcf ${params.faf} ${params.af_cutoff}
    bcftools view slivar_vcf.vcf -Oz -o ${input_vcf.SimpleName}.vcf.gz
    """
}

workflow preanalysis {

  // Pre-analysis subworkflow

  take:
    input_vcfs
    chr_format
    fasta_file
    gnomad_annotation
    normalization_script
    freq_filter_script
 
  main: 
    input_files = format_vcf_files(input_vcfs.combine(chr_format))
    fasta_refs = format_reference_files(fasta_file)
    fasta_ref = fasta_refs.fa
    normalized_files = normalize_files(input_files, fasta_ref, normalization_script)

    if ( params.faf ) {

        formatted_vcfs = slivar_annotation(normalized_files.combine(gnomad_annotation).combine(fasta_ref).combine(freq_filter_script))

    } else {

        formatted_vcfs = normalized_files

    }

  emit: 
    formatted_vcfs
    fasta_ref

}