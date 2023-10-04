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
    
    stageInMode 'copy'
    label 'inParallel'

    input:
        path fasta_file

    output:
        path "${fasta_file.SimpleName}.fa", emit: fa
        path "${fasta_file.SimpleName}.fa.fai", emit: fai

    script:
    """
    if file -b --mime-type ${fasta_file} | grep -q gzip; then
        gunzip -c ${fasta_file.SimpleName}.fa.gz > ${fasta_file.SimpleName}.fa
    else
        true
    fi
    sed -i 's/>chr/>/g' ${fasta_file.SimpleName}.fa > ${fasta_file.SimpleName}.fa
    samtools faidx ${fasta_file.SimpleName}.fa
    """

}

process normalize_files {

    /* Function to normalize and index vcfs using bcftools */ 
    
    input:
        each vcf_file
        path fasta_file

    output:
        tuple path("${vcf_file.SimpleName}.vcf.gz"), path("${vcf_file.SimpleName}.vcf.gz.tbi")

    script:
    """
    bcftools norm -cs -f ${fasta_file} ${vcf_file} -Oz -o normalized_vcf.vcf.gz
    mv normalized_vcf.vcf.gz ${vcf_file.SimpleName}.vcf.gz
    bcftools index -f -t ${vcf_file.SimpleName}.vcf.gz
    """

}

process annotate_allele_frequency {

    /* Uses Slivar to annotate Gnomad allele
    frequencies */
    
    label 'inParallel'

    input:
        tuple path(input_vcf), path(input_tbi), path(gnomad_annotation), path(fasta_bgzip), path(fasta_fai)

    output:
        tuple path("${input_vcf.SimpleName}.vcf.gz"), path("${input_tbi.SimpleName}.vcf.gz.tbi")

    script:
    """
    slivar expr -g ${gnomad_annotation} --vcf ${input_vcf} > ${input_vcf}
    """
}

workflow preanalysis {

  // Pre-analysis subworkflow

  take:
    input_vcfs
    chr_format
    fasta_file
    gnomad_annotation
 
  main: 
    input_files = format_vcf_files(input_vcfs.combine(chr_format))
    fasta_refs = format_reference_files(fasta_file)
    fasta_ref = fasta_refs.fa
    normalized_files = normalize_files(input_files, fasta_refs.fa)

    if ( params.faf ) {

        formatted_vcfs = annotate_allele_frequency(normalized_files.combine(gnomad_annotation).combine(fasta_refs.fa).combine(fasta_refs.fai))

    } else {

        formatted_vcfs = normalized_files

    }

  emit: 
    formatted_vcfs
    fasta_ref

}