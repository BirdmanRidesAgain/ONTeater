#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
/*
 * Keiler Collier
 * ONTeater V1 - used to run genome assembly from concatenated, Kraken-filtered readfiles.
 *
 * Started: 21 Feb 2024
 * Last update: 12 Nov 2025
 */


//include { REMOVE_CONTAMINANTS } from './modules/remove_contaminants_kraken2/main.nf' //initial filtering for contaminants
// Import modules
include { PRINT_HELP } from './modules/print_help/main.nf'
include { PREPROCESS_DATA } from './modules/preprocess_reads.nf'
include { PRIMARY_ASSEMBLY } from './modules/primary_assembly.nf'
include { POSTPROCESS_ASSEMBLY } from './modules/postprocess_assembly.nf'
include { QC_WORKFLOW } from './modules/qc_workflow.nf'

workflow {
    main:
    def ont_reads = params.ONT_rds
    def workflow_mode = params.workflow ?: 'run'
    def valid_modes = ['run', 'trim', 'assemble', 'postprocess', 'qc']

    if (params.help) {
        PRINT_HELP()
        exit 10
    }

    if (!valid_modes.contains(workflow_mode)) {
        error "Unsupported --workflow '${workflow_mode}'. Valid options: ${valid_modes.join(', ')}"
    }

    if (['run', 'trim', 'assemble'].contains(workflow_mode) && !ont_reads) {
        error "No ONT reads provided. Use --ONT_rds."
    }
    if (workflow_mode == 'postprocess' && (!params.flye_asm || !params.nd_asm)) {
        error "Workflow mode 'postprocess' requires --flye_asm and --nd_asm."
    }
    if (workflow_mode == 'qc' && !params.final_asm) {
        error "Workflow mode 'qc' requires --final_asm."
    }

    log.info """\
        O N T E A T E R - N F   P I P E L I N E
        ===================================
        Project directory       : $projectDir
        Input ONT longreads     : ${ont_reads ?: 'N/A'}
        Workflow mode           : ${workflow_mode}
        Genome size             : ${params.genome_size}
        BUSCO lineage           : ${params.BUSCO_lineage}
        """
        .stripIndent()

    def ch_rawreads = ont_reads ? channel.fromPath(ont_reads) : channel.empty()

    def ch_trimreads
    def ch_raw_viz
    def ch_trim_viz
    if (['run', 'trim'].contains(workflow_mode)) {
        PREPROCESS_DATA(ch_rawreads)
        ch_trimreads = PREPROCESS_DATA.out.trimreads
        ch_raw_viz = PREPROCESS_DATA.out.rawread_viz
        ch_trim_viz = PREPROCESS_DATA.out.trimread_viz
    } else {
        // assemble mode starts directly from provided reads.
        ch_trimreads = ch_rawreads
        ch_raw_viz = channel.empty()
        ch_trim_viz = channel.empty()
    }

    def ch_trimreads_publish = (['run', 'trim'].contains(workflow_mode)) ? ch_trimreads : channel.empty()
    def ch_flye_polished = channel.empty()
    def ch_nextdenovo_polished = channel.empty()
    def ch_merged_assembly = channel.empty()
    def ch_purged_assembly = channel.empty()
    def ch_qc_quast = channel.empty()
    def ch_qc_compleasm = channel.empty()

    if (['run', 'assemble'].contains(workflow_mode)) {
        PRIMARY_ASSEMBLY(ch_trimreads)
        ch_flye_polished = PRIMARY_ASSEMBLY.out.flye_polished
        ch_nextdenovo_polished = PRIMARY_ASSEMBLY.out.nextdenovo_polished
    }

    if (workflow_mode == 'postprocess') {
        def flye_path = file(params.flye_asm)
        def nd_path = file(params.nd_asm)
        def flye_sample = flye_path.simpleName.replaceFirst(/_Flye(_racon)?$/, '')
        def nd_sample = nd_path.simpleName.replaceFirst(/_nextDenovo(_racon)?$/, '')
        ch_flye_polished = channel.of(tuple(flye_sample, 'Flye', flye_path))
        ch_nextdenovo_polished = channel.of(tuple(nd_sample, 'nextDenovo', nd_path))
    }

    if (['run', 'postprocess'].contains(workflow_mode)) {
        POSTPROCESS_ASSEMBLY(ch_flye_polished, ch_nextdenovo_polished)
        ch_merged_assembly = POSTPROCESS_ASSEMBLY.out.merged_assembly
        ch_purged_assembly = POSTPROCESS_ASSEMBLY.out.purged_assembly
    }

    if (['run', 'qc'].contains(workflow_mode)) {
        def ch_qc_input = (workflow_mode == 'qc')
            ? channel.of(tuple(file(params.final_asm).simpleName, file(params.final_asm)))
            : ch_purged_assembly.map { sample_id, fasta -> tuple(sample_id, fasta) }
        QC_WORKFLOW(ch_qc_input)
        ch_qc_quast = QC_WORKFLOW.out.quast
        ch_qc_compleasm = QC_WORKFLOW.out.compleasm
    }

    publish:
    raw_read_qc = ch_raw_viz
    trim_read_qc = ch_trim_viz
    trimmed_reads = ch_trimreads_publish
    flye_polished = ch_flye_polished
    nextdenovo_polished = ch_nextdenovo_polished
    merged_assembly = ch_merged_assembly
    purged_assembly = ch_purged_assembly
    qc_quast = ch_qc_quast
    qc_compleasm = ch_qc_compleasm
}

output {
    raw_read_qc {
        path "${params.prefix}/reads"
    }
    trim_read_qc {
        path "${params.prefix}/reads"
    }
    trimmed_reads {
        path "${params.prefix}/reads"
    }
    flye_polished {
        path "${params.prefix}/assemblies"
    }
    nextdenovo_polished {
        path "${params.prefix}/assemblies"
    }
    merged_assembly {
        path "${params.prefix}/assemblies"
    }
    purged_assembly {
        path "${params.prefix}/assemblies"
    }
    qc_quast {
        path "${params.prefix}/assemblies"
    }
    qc_compleasm {
        path "${params.prefix}/assemblies"
    }
}
