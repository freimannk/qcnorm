#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process normalise_microarray{
    publishDir "${params.outdir}/$study_id/normalised/microarray", mode: 'copy', pattern: "qtl_group_split_norm/*"
    
    container = 'quay.io/eqtlcatalogue/eqtlutils:v22.11.1'
    
    input:
    tuple val(study_id), val(sample_group), file(quant_results_path), file(sample_metadata)
    path pheno_metadata
    
    output:
    tuple val(study_id), file("qtl_group_split_norm/*"), emit: qtlmap_tsv_input_ch

    script:
    filter_qc = params.norm_filter_qc ? "--filter_qc TRUE" : ""
    keep_XY = params.norm_keep_XY ? "--keep_XY TRUE" : ""
    eqtl_utils_path = params.eqtl_utils_path ? "--eqtlutils ${params.eqtl_utils_path}" : ""
    """
    Rscript $baseDir/bin/normalisation/normaliseCountMatrix.R\
      -c $quant_results_path\
      -s $sample_metadata\
      -p $pheno_metadata\
      -n $study_id\
      -o .\
      -q HumanHT-12_V4\
      $filter_qc\
      $keep_XY\
      $eqtl_utils_path

    """
}

process normalise_RNAseq_ge{
    publishDir "${params.outdir}/$study_id/normalised/ge", mode: 'copy', pattern: "norm_not_filtered/*"
    publishDir "${params.outdir}/$study_id/normalised/ge", mode: 'copy', pattern: "qtl_group_split_norm/*"
    publishDir "${params.outdir}/$study_id/normalised/ge", mode: 'copy', pattern: "qtl_group_split_norm_anonym/*"
    publishDir "${params.outdir}/$study_id/normalised/ge", mode: 'copy', pattern: "per_million_normalised/*"
    publishDir "${params.outdir}/$study_id/normalised/ge", mode: 'copy', pattern: "qtl_group_median_tpms/*"
    publishDir "${params.outdir}/$study_id/normalised/", mode: 'copy', pattern: "*_tpm.tsv.gz"

    container = 'quay.io/eqtlcatalogue/eqtlutils:v22.11.1'
    
    input:
    tuple val(study_id), val(sample_group), file(quant_results_path), file(sample_metadata)
    path pheno_metadata
    
    output:
    path "norm_not_filtered/*"
    path "*_95quantile_tpm.tsv.gz", emit: quantile_tpm_file
    tuple val(study_id),val(sample_group), file("*_median_tpm.tsv.gz"), emit: median_tpm_file
    tuple val(study_id),val(sample_group), file("qtl_group_split_norm/*"), emit: qtlmap_tsv_input_ch
    tuple val(study_id), val(sample_group), file(quant_results_path), file(sample_metadata), file("*_95quantile_tpm.tsv.gz"), emit: inputs_with_quant_tpm_ch
    path "qtl_group_split_norm_anonym/*"
    path "per_million_normalised/*"
    path "qtl_group_median_tpms/*"

    script:
    filter_qc = params.norm_filter_qc ? "--filter_qc TRUE" : ""
    keep_XY = params.norm_keep_XY ? "--keep_XY TRUE" : ""
    eqtl_utils_path = params.eqtl_utils_path ? "--eqtlutils ${params.eqtl_utils_path}" : ""
    """
    Rscript $baseDir/bin/normalisation/normaliseCountMatrix.R\
      -c $quant_results_path/featureCounts/merged_gene_counts.tsv.gz\
      -s $sample_metadata\
      -p $pheno_metadata\
      -n $study_id\
      -o .\
      -q gene_counts\
      $filter_qc\
      $keep_XY\
      $eqtl_utils_path

    """
}
process format_majiq_inputs{
    publishDir "${params.outdir}/$study_id/normalised/majiq", mode: 'copy', pattern: '*majiq_metadata.tsv.gz'

    container = 'quay.io/kfkf33/polars_qcnorm:v1.0'
    
    input:
    tuple val(study_id), val(sample_group), file(quant_results_path), file(sample_metadata)
    path pheno_metadata
    
    output:
    tuple val(study_id),val(sample_group), file("${study_id}.${sample_group}_pheno_matrix.tsv.gz"), emit: phenotype_matrix
    tuple val(study_id),val(sample_group), file("${study_id}.${sample_group}_majiq_metadata.tsv.gz"), emit: phenotype_metadata
    
    script:
    filter_qc = params.norm_filter_qc ? "-q TRUE" : ""
    keep_XY = params.norm_keep_XY ? "-x TRUE" : ""
    """
    $baseDir/bin/normalisation/format_majiq_inputs.py \
        -p $pheno_metadata \
        -m $quant_results_path/majiq/majiq_quantified_psis/majiq_quantified_psis.tsv.gz \
        -s $sample_metadata \
        -o ${study_id}.${sample_group} \
        $keep_XY \
        $filter_qc

    """
}

process normalise_RNAseq_majiq{
    publishDir "${params.outdir}/$study_id/normalised/majiq", mode: 'copy'
    container = "quay.io/eqtlcatalogue/eqtlutils:v25.10.1"
    
    input:
    tuple val(study_id), val(sample_group), file(phenotype_matrix)
    
    output:
     tuple val(study_id),val(sample_group), file("${study_id}.${sample_group}.tsv.gz"), emit: qtlmap_tsv_input_ch

    script:
      min_experiments = params.majiq_min_experiments ? "-m ${params.majiq_min_experiments}" : ""
   
    """
    Rscript $baseDir/bin/normalisation/normaliseMajiq.R\
      -p $phenotype_matrix\
      -o ${study_id}.${sample_group}\
      $min_experiments


    """
}

process normalise_RNAseq_exon{
    publishDir "${params.outdir}/$study_id/normalised/exon", mode: 'copy'
    
    container = 'quay.io/eqtlcatalogue/eqtlutils:v22.11.1'
    
    input:
    tuple val(study_id), val(sample_group), file(quant_results_path), file(sample_metadata), file(tpm_quantile)
    path pheno_metadata

    output:
    path "norm_not_filtered/*"
    tuple val(study_id),val(sample_group), file("qtl_group_split_norm/*"), emit: qtlmap_tsv_input_ch
    path "per_million_normalised/*"

    script:
    filter_qc = params.norm_filter_qc ? "--filter_qc TRUE" : ""
    keep_XY = params.norm_keep_XY ? "--keep_XY TRUE" : ""
    eqtl_utils_path = params.eqtl_utils_path ? "--eqtlutils ${params.eqtl_utils_path}" : ""
    """
    Rscript $baseDir/bin/normalisation/normaliseCountMatrix.R\
      -c $quant_results_path/dexseq_exon_counts/merged_exon_counts.tsv.gz\
      -s $sample_metadata\
      -p $pheno_metadata\
      -n $study_id\
      -o .\
      -q exon_counts\
      -t $tpm_quantile\
      $filter_qc\
      $keep_XY\
      $eqtl_utils_path

    """
}

process normalise_RNAseq_tx{
    publishDir "${params.outdir}/$study_id/normalised/tx", mode: 'copy'
    
    container = 'quay.io/eqtlcatalogue/eqtlutils:v22.11.1'
    
    input:
    tuple val(study_id), val(sample_group), file(quant_results_path), file(sample_metadata), file(tpm_quantile)
    path pheno_metadata
    
    output:
    path "norm_not_filtered/*"
    tuple val(study_id), val(sample_group), file("qtl_group_split_norm/*"), emit: qtlmap_tsv_input_ch
    path "per_million_normalised/*"

    script:
    filter_qc = params.norm_filter_qc ? "--filter_qc TRUE" : ""
    keep_XY = params.norm_keep_XY ? "--keep_XY TRUE" : ""
    eqtl_utils_path = params.eqtl_utils_path ? "--eqtlutils ${params.eqtl_utils_path}" : ""
    """
    Rscript $baseDir/bin/normalisation/normaliseCountMatrix.R\
      -c $quant_results_path/Salmon/merged_counts/TPM/gencode.v39.transcripts.TPM.merged.tsv.gz\
      -s $sample_metadata\
      -p $pheno_metadata\
      -n $study_id\
      -o .\
      -q transcript_usage\
      -t $tpm_quantile\
      $filter_qc\
      $keep_XY\
      $eqtl_utils_path

    """
}

process normalise_RNAseq_txrev{
    publishDir "${params.outdir}/$study_id/normalised/txrev", mode: 'copy'
    
    container = 'quay.io/eqtlcatalogue/eqtlutils:v22.11.1'
    
    input:
    tuple val(study_id), val(sample_group), file(quant_results_path), file(sample_metadata), file(tpm_quantile)
    path pheno_metadata
    
    output:
    path "norm_not_filtered/*"
    tuple val(study_id),val(sample_group), file("qtl_group_split_norm/*"), emit: qtlmap_tsv_input_ch
    path "per_million_normalised/*"

    script:
    filter_qc = params.norm_filter_qc ? "--filter_qc TRUE" : ""
    keep_XY = params.norm_keep_XY ? "--keep_XY TRUE" : ""
    eqtl_utils_path = params.eqtl_utils_path ? "--eqtlutils ${params.eqtl_utils_path}" : ""
    """
    Rscript $baseDir/bin/normalisation/normaliseCountMatrix.R\
      -c $quant_results_path/Salmon/merged_counts/TPM/\
      -s $sample_metadata\
      -p $pheno_metadata\
      -n $study_id\
      -o .\
      -q txrevise\
      -t $tpm_quantile\
      $filter_qc\
      $keep_XY\
      $eqtl_utils_path

    """
}


process normalise_RNAseq_leafcutter{
    publishDir "${params.outdir}/$study_id/normalised/leafcutter", mode: 'copy'
    
    label 'process_medium'
    container = 'quay.io/eqtlcatalogue/eqtlutils:v22.11.1'
    
    input:
    tuple val(study_id), val(sample_group), file(quant_results_path), file(sample_metadata), file(tpm_quantile)
    path transcript_meta
    path intron_annotation
    
    output:
    path "norm_not_filtered/*"
    tuple val(study_id),val(sample_group), file("qtl_group_split_norm/*"), emit: qtlmap_tsv_input_ch
    path "${sample_group}_leafcutter_metadata.txt.gz", emit: leafcutter_metadata
    path "per_million_normalised/*"

    script:
    filter_qc = params.norm_filter_qc ? "--filter_qc TRUE" : ""
    keep_XY = params.norm_keep_XY ? "--keep_XY TRUE" : ""
    eqtl_utils_path = params.eqtl_utils_path ? "--eqtlutils ${params.eqtl_utils_path}" : ""
    """
    # Make leafcutter phenotype metadata file
    Rscript $baseDir/bin/normalisation/makeLeafcutterMetadata.R\
      -c $quant_results_path/leafcutter/leafcutter_perind_numers.counts.formatted.gz\
      -t $transcript_meta\
      -i $intron_annotation\
      -o ${sample_group}_leafcutter_metadata.txt.gz\
      $eqtl_utils_path

    Rscript $baseDir/bin/normalisation/normaliseCountMatrix.R\
      -c $quant_results_path/leafcutter/leafcutter_perind_numers.counts.formatted.gz\
      -s $sample_metadata\
      -p ${sample_group}_leafcutter_metadata.txt.gz\
      -n $study_id\
      -o .\
      -q leafcutter\
      -t $tpm_quantile\
      $filter_qc\
      $keep_XY\
      $eqtl_utils_path

    """
}