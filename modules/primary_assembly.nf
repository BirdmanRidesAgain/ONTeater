include { ASSEMBLE_FLYE } from './assemble_flye/main.nf'
include { GET_NEXTDENOVO_PARAMS } from './get_nextdenovo_params/main.nf'
include { ASSEMBLE_NEXTDENOVO } from './assemble_nextdenovo/main.nf'
include { POLISH_RACON as POLISH_RACON_FLYE } from './polish_racon/main.nf'
include { POLISH_RACON as POLISH_RACON_ND } from './polish_racon/main.nf'

workflow PRIMARY_ASSEMBLY {

    take:
    ch_trimreads

    main:
    // Convert trimmed read paths into a consistent metadata tuple expected by assembly modules.
    ch_reads_meta = ch_trimreads.map { reads ->
        def sample_id = reads.simpleName.replaceFirst(/_trim$/, '')
        tuple(sample_id, 'trim', reads)
    }

    ASSEMBLE_FLYE(ch_reads_meta)
    ch_flye = ASSEMBLE_FLYE.out

    // Build nextDenovo config from template and current genome-size setting.
    ch_nd_template = channel.value(file("${projectDir}/bin/nd_run.cfg"))
    GET_NEXTDENOVO_PARAMS(ch_flye, ch_nd_template)
    ch_nd_cfg = GET_NEXTDENOVO_PARAMS.out

    ASSEMBLE_NEXTDENOVO(ch_reads_meta, ch_nd_cfg)
    ch_nextdenovo = ASSEMBLE_NEXTDENOVO.out

    ch_polish_flye_input = ch_flye
        .join(ch_reads_meta)
        .map { sample_id, assembler, fasta, trim_status, reads -> tuple(sample_id, assembler, fasta, trim_status, reads) }

    ch_polish_nd_input = ch_nextdenovo
        .join(ch_reads_meta)
        .map { sample_id, assembler, fasta, trim_status, reads -> tuple(sample_id, assembler, fasta, trim_status, reads) }

    POLISH_RACON_FLYE(ch_polish_flye_input)
    POLISH_RACON_ND(ch_polish_nd_input)

    emit:
    flye_polished = POLISH_RACON_FLYE.out
    nextdenovo_polished = POLISH_RACON_ND.out
    polished_assemblies = POLISH_RACON_FLYE.out.mix(POLISH_RACON_ND.out)
}
