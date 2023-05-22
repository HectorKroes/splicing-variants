#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SpliceAI GPU subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Implements the usage of SpliceAI v1.3.1 with the option 
    of using precalculated scores to cut processing time
    and save computing resources
----------------------------------------------------------------------------------------
*/

process annotate_precalculated_scores {

  /* Function receives input vcf, and precalculated scores,
  annotates all the variants for which there are available
  scores and then divides file in two: pcs (which contains
  precalculated scores) and tbc (which contains variants
  not annotated, that will be redirected to receive de novo
  predictions by SpliceAI) */
  
  label 'inParallel'

  input:
    tuple path(input_vcf), path(index_file)
    path indel_annotation
    path indel_annotation_index
    path snv_annotation
    path snv_annotation_index

  output:
    tuple val("${input_vcf.baseName}"), path("pcs_${input_vcf.baseName}"), emit: pcs_ch
    tuple val("${input_vcf.baseName}"), path("tbc_${input_vcf.baseName}"), emit: tbc_ch

  script:
    """
    bcftools annotate -c 'INFO' -a ${indel_annotation} ${input_vcf} -O z -o pcs1.vcf.gz
    bcftools index -f -t pcs1.vcf.gz
    bcftools annotate -c 'INFO' -a ${snv_annotation} pcs1.vcf.gz -O z -o pcs2.vcf.gz
    zcat pcs2.vcf.gz | grep '^#' > pcs_${input_vcf.baseName}
    zcat pcs2.vcf.gz | grep SpliceAI= >> pcs_${input_vcf.baseName}
    zcat pcs2.vcf.gz | grep -v SpliceAI= > tbc_${input_vcf.baseName}
    """
}

process predict_de_novo_variants {

  /* Function receives tuple containing the tbc file, 
  the equivalent file basename and a fasta reference
  for running SpliceAI, outputting the results as a tuple
  of the file basename and a dnv file containing the de 
  novo predictions. The function also ensures that the fasta
  file has compatible chr formats with the files presented. */
  
  stageInMode 'copy'
  label 'inSeries'
  
  input:
    tuple val(basename), path(to_be_computed_vcf), path(ref_fasta)

  output:
    tuple val("${basename}"), path("dnv_${basename}")

  script:
    """
    spliceai -I tbc_${basename} -O dnv_${basename} -R ${ref_fasta} -A ${params.ref} -D ${params.spliceai_max_length}
    """
}

process fuse_temporary_vcfs {
  
  /* Function that gets dnv and pcs vcf files and concatenates
  them to create a result file with all variants */

  label 'inParallel'
  
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

process no_pcv_adequation {

  /* Function that fits files that will not be annotated 
  by precalculated scores into a tuple compatible with 
  predict_de_novo_variants function input format */ 

  label 'inParallel'

  input:
    tuple path(input_vcf), path(index_file)

  output:
    tuple val("${input_vcf.baseName}"), path("tbc_${input_vcf.baseName}")

  script:
    """
    gunzip -c ${input_vcf} > tbc_${input_vcf.baseName}
    """

}

workflow spliceai_gpu {

  /* SpliceAI GPU subworkflow - pcv parameter check determines 
  whether to use precalculated scores or not */

  take:
    input_files
    snv_annotation
    indel_annotation
    snv_annotation_index
    indel_annotation_index
    fasta_ref
 
  main: 

    if ( params.pcv ) {

      pcs_channel = annotate_precalculated_scores(input_files, indel_annotation, indel_annotation_index, snv_annotation, snv_annotation_index)

      dnv_predictions = predict_de_novo_variants(pcs_channel.tbc_ch.combine(fasta_ref))

      temporary_vcfs = (pcs_channel.pcs_ch).join(dnv_predictions)

      spliceai_results = fuse_temporary_vcfs(temporary_vcfs)

    } else {

      tbc_channel = no_pcv_adequation(input_files)

      spliceai_results = predict_de_novo_variants(tbc_channel.combine(fasta_ref))
    
    }

  emit: 
    spliceai_results
  
}