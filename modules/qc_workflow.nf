include { QC_QUAST } from './qc_processes.nf'
include { QC_COMPLEASM } from './qc_processes.nf'

workflow QC_WORKFLOW {
    take:
    ch_assembly

    main:
    QC_QUAST(ch_assembly)
    QC_COMPLEASM(ch_assembly)

    emit:
    quast = QC_QUAST.out
    compleasm = QC_COMPLEASM.out
}
