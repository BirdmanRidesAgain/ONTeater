include { VISUALIZE_READS_NANOPLOT as VIZ_NANO_RAW } from './visualize_reads_nanoplot'
include { VISUALIZE_READS_NANOPLOT as VIZ_NANO_TRIM } from './visualize_reads_nanoplot'
include { TRIM_READS_CHOPPER } from './trim_reads_chopper' 

workflow PREPROCESS_READS {

    take:
    ch_reads

    main:
    VIZ_NANO_RAW(ch_reads, 'raw')
    TRIM_READS_CHOPPER(ch_reads)
    VIZ_NANO_TRIM(TRIM_READS_CHOPPER.out.trimreads, 'filtered')

    emit:
    trimreads = TRIM_READS_CHOPPER.out.trimreads
    rawread_viz = VIZ_NANO_RAW.out.viz
    trimread_viz = VIZ_NANO_TRIM.out.viz
}

workflow PREPROCESS_DATA {
    take:
    ch_reads

    main:
    PREPROCESS_READS(ch_reads)

    emit:
    trimreads = PREPROCESS_READS.out.trimreads
    rawread_viz = PREPROCESS_READS.out.rawread_viz
    trimread_viz = PREPROCESS_READS.out.trimread_viz
}