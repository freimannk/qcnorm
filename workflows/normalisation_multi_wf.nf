#!/usr/bin/env nextflow
nextflow.enable.dsl=2

qtl_inputs_def_file = file("$baseDir/assets/qtlmap_inputs.tsv")
qtl_inputs_def_file.copyTo("${params.outdir}/qtlmap_inputs.tsv")
qtlmap_inputs_file = file("${params.outdir}/qtlmap_inputs.tsv")
pheno_metadata_list = [
    "ge": params.ge_pheno_meta_path,
    "exon": params.exon_pheno_meta_path,
    "tx": params.tx_pheno_meta_path,
    "txrev": params.txrev_pheno_meta_path,
    "microarray": params.array_pheno_meta_path
]

dataset_map = file(params.dataset_info)
    .splitCsv(header:true, sep:'\t')
    .collectEntries { row ->
        [ "${row.study_id.toString()}|${row.quant_method.toString()}|${row.sample_group.toString()}" : row.dataset_id.toString() ]
    }
    .collectEntries { k,v -> [ k.toString(), v.toString() ] }

include { normalise_microarray; normalise_RNAseq_ge ; normalise_RNAseq_exon ; normalise_RNAseq_tx ; normalise_RNAseq_txrev; normalise_RNAseq_leafcutter; format_majiq_inputs; normalise_RNAseq_majiq } from  '../modules/normalisation_multi'

def add_to_qtlmap_input_tsv(qtlgroup_quantiletpm_ch, quant_method) {
    def method_paths = [
        microarray  : { item -> [normalized_matrix: "${params.outdir}/${item[0]}/normalised/${quant_method}/qtl_group_split_norm/${item[4].fileName}", phenotype_metadata: 'null.txt'] },
        leafcutter  : { item -> [
                            normalized_matrix: "${params.outdir}/${item[0]}/normalised/${quant_method}/qtl_group_split_norm/${item[5].fileName}",
                            phenotype_metadata: "${params.outdir}/${item[0]}/normalised/${quant_method}/${item[1]}_leafcutter_metadata.txt.gz",
                            median_tpm: "${params.outdir}/${item[0]}/normalised/ge/qtl_group_median_tpms/${item[0]}_ge_${item[1]}_median_tpm.tsv.gz"
                        ] },
        majiq       : { item -> [
                            normalized_matrix: "${params.outdir}/${item[0]}/normalised/${quant_method}/qtl_group_split_norm/${item[0]}.${item[1]}.tsv.gz",
                            phenotype_metadata: "${params.outdir}/${item[0]}/normalised/${quant_method}/${item[0]}.${item[1]}_majiq_metadata.tsv.gz",
                            median_tpm: "${params.outdir}/${item[0]}/normalised/ge/qtl_group_median_tpms/${item[0]}_ge_${item[1]}_median_tpm.tsv.gz"
                        ] },
        default     : { item -> [
                            normalized_matrix: "${params.outdir}/${item[0]}/normalised/${quant_method}/qtl_group_split_norm/${item[5].fileName}",
                            phenotype_metadata: pheno_metadata_list[quant_method],
                            median_tpm: "${params.outdir}/${item[0]}/normalised/ge/qtl_group_median_tpms/${item[0]}_ge_${item[1]}_median_tpm.tsv.gz"
                        ] }
    ]

    def paths_fn = method_paths.get(quant_method, method_paths.default)

    qtlgroup_quantiletpm_ch
        .map { item ->
            def key = "${item[0]}|${quant_method}|${item[1]}"
            def dataset_id = dataset_map[key]
            def paths = paths_fn(item)
            tuple(item, dataset_id, paths)
        }
        .collectFile(storeDir: "${params.outdir}/qtl_group_inputs") { x ->
            def item = x[0]
            def dataset_id = x[1]
            def paths = x[2]

            [
                "${item[0]}_${item[1]}_${quant_method}_tsv_inputs.txt",
                "${dataset_id}\t" +
                "${paths.normalized_matrix}\t" +
                (paths.phenotype_metadata ?: '') + "\t" +
                "${item[2]}\t" +
                "${item[3]}\t" +
                (paths.median_tpm ?: 'null.txt') + '\n'
            ]
        }
        .subscribe { qtlmap_inputs_file.append(it.text) }
}

workflow {
    normalise()
}

workflow normalise {
        //study_id	sample_group	quant_results_path	sample_meta_path	vcf_file
    Channel.fromPath(params.input_tsv)
        .ifEmpty { error "Cannot find input_tsv file in: ${params.input_tsv}" }
        .splitCsv(header: true, sep: '\t', strip: true)
        .map{row -> [ row.study_id, row.sample_group, file(row.quant_results_path), file(row.sample_meta_path) ]}
        .set { study_file_ge_ch }

    Channel.fromPath(params.input_tsv)
        .ifEmpty { error "Cannot find input_tsv file in: ${params.input_tsv}" }
        .splitCsv(header: true, sep: '\t', strip: true)
        .map{row -> [ row.study_id, row.sample_group, row.sample_meta_path, row.vcf_file ]}
        .set { output_tsv_ch }
        
    if (params.is_microarray){
        normalise_microarray(study_file_ge_ch, Channel.fromPath(params.array_pheno_meta_path, checkIfExists: true).collect())

        add_to_qtlmap_input_tsv(output_tsv_ch
            .combine(normalise_microarray.out.qtlmap_tsv_input_ch, by: 0).transpose(), "microarray")
    }
    else {
        normalise_RNAseq_ge(study_file_ge_ch, Channel.fromPath(params.ge_pheno_meta_path, checkIfExists: true).collect())
        
        add_to_qtlmap_input_tsv(output_tsv_ch
            .join(normalise_RNAseq_ge.out.median_tpm_file, by: [0,1])
            .combine(normalise_RNAseq_ge.out.qtlmap_tsv_input_ch, by: [0,1]).transpose(), "ge")
        
        if (!params.skip_exon_norm) {
            normalise_RNAseq_exon(normalise_RNAseq_ge.out.inputs_with_quant_tpm_ch, Channel.fromPath(params.exon_pheno_meta_path, checkIfExists: true).collect())
        
            add_to_qtlmap_input_tsv(output_tsv_ch
                .join(normalise_RNAseq_ge.out.median_tpm_file, by: [0,1])
                .combine(normalise_RNAseq_exon.out.qtlmap_tsv_input_ch, by: [0,1]).transpose(), "exon")
        }
        if (!params.skip_tx_norm) {
            normalise_RNAseq_tx(normalise_RNAseq_ge.out.inputs_with_quant_tpm_ch, Channel.fromPath(params.tx_pheno_meta_path, checkIfExists: true).collect())

            add_to_qtlmap_input_tsv(output_tsv_ch
                .join(normalise_RNAseq_ge.out.median_tpm_file, by: [0,1])
                .combine(normalise_RNAseq_tx.out.qtlmap_tsv_input_ch, by: [0,1]).transpose(), "tx")
        } 
        if (!params.skip_txrev_norm) {
            normalise_RNAseq_txrev(normalise_RNAseq_ge.out.inputs_with_quant_tpm_ch, Channel.fromPath(params.txrev_pheno_meta_path, checkIfExists: true).collect())

            add_to_qtlmap_input_tsv(output_tsv_ch
                .join(normalise_RNAseq_ge.out.median_tpm_file, by: [0,1])
                .combine(normalise_RNAseq_txrev.out.qtlmap_tsv_input_ch, by: [0,1]).transpose(), "txrev")
        }
        if (!params.skip_leafcutter_norm) {
            normalise_RNAseq_leafcutter(normalise_RNAseq_ge.out.inputs_with_quant_tpm_ch, 
                                        Channel.fromPath(params.leafcutter_transcript_meta, checkIfExists: true).collect(),
                                        Channel.fromPath(params.leafcutter_intron_annotation, checkIfExists: true).collect())


            add_to_qtlmap_input_tsv(output_tsv_ch
                .join(normalise_RNAseq_ge.out.median_tpm_file, by: [0,1])
                .combine(normalise_RNAseq_leafcutter.out.qtlmap_tsv_input_ch, by: [0,1]).transpose(), "leafcutter")
        }

        if (!params.skip_majiq_norm) {
            format_majiq_inputs(study_file_ge_ch, Channel.fromPath(params.ge_pheno_meta_path, checkIfExists: true).collect())
            normalise_RNAseq_majiq(format_majiq_inputs.out.phenotype_matrix)

            add_to_qtlmap_input_tsv(output_tsv_ch
                .join(normalise_RNAseq_ge.out.median_tpm_file, by: [0,1])
                .combine(normalise_RNAseq_majiq.out.qtlmap_tsv_input_ch, by: [0,1]).transpose(), "majiq")
        }
    }
}
