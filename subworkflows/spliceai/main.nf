#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SpliceAI subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Implements the usage of SpliceAI v1.3.1 with the option 
    of using precalculated scores to cut processing time
----------------------------------------------------------------------------------------
*/

process spliceai_annotate_precalculated_scores {

  /* Function receives input vcf, and precalculated scores,
  annotates all the variants for which there are available
  scores and then divides file in two: pcs (which contains
  precalculated scores) and tbc (which contains variants
  not annotated, that will be redirected to receive de novo
  predictions by SpliceAI) */
  
  tag "${input_vcf.baseName}"

  input:
    each input_vcf
    path indel_annotation
    path indel_annotation_index
    path snv_annotation
    path snv_annotation_index

  output:
    tuple val("${input_vcf.baseName}"), path("pcs_${input_vcf.baseName}"), emit: pcs_ch
    tuple val("${input_vcf.baseName}"), path("tbc_${input_vcf.baseName}"), emit: tbc_ch

  script:
    """
    bcftools index -f -t ${input_vcf}
    bcftools annotate -c 'INFO' -a ${indel_annotation} ${input_vcf} -O z -o pcs1.vcf.gz
    bcftools index -f -t pcs1.vcf.gz
    bcftools annotate -c 'INFO' -a ${snv_annotation} pcs1.vcf.gz -O z -o pcs2.vcf.gz
    zcat pcs2.vcf.gz | grep '^#' > pcs_${input_vcf.baseName}
    zcat pcs2.vcf.gz | grep SpliceAI= >> pcs_${input_vcf.baseName}
    zcat pcs2.vcf.gz | grep -v SpliceAI= > tbc_${input_vcf.baseName}
    """
}

process spliceai_predict_de_novo_variants {

  /* Function receives tuple containing the tbc file, 
  the equivalent file basename and a fasta reference
  for running SpliceAI, outputting the results as a tuple
  of the file basename and a dnv file containing the de 
  novo predictions. */
  
  tag "${basename}"
  cpus params.t
  
  input:
    tuple val(basename), path(to_be_computed_vcf), path(ref_fasta)

  output:
    tuple val("${basename}"), path("dnv_${basename}")

  script:
    """
    spliceai -I tbc_${basename} -O dnv_${basename} -R ${ref_fasta} -A grch38
    """
}

process spliceai_fuse_temporary_vcfs {
  
  /* Function get dnv and pcs vcf files and concatenates
  them to create a result file with all variants */

  tag "${basename}"
  
  input:
    tuple val(basename), path(vcf1), path(vcf2)

  output:
    path "${basename}.gz"

  script:
    """
    bcftools view ${vcf1} -Oz -o ${vcf1}.gz
    bcftools index -f -t ${vcf1}.gz
    bcftools view ${vcf2} -Oz -o ${vcf2}.gz
    bcftools index -f -t ${vcf2}.gz
    bcftools concat -a -O z -o ${basename}.gz ${vcf1}.gz ${vcf2}.gz
    """

}