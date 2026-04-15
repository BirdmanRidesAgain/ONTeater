include { ASSESS_ASSEMBLY_QUAST } from './assess_assembly_quast/main.nf'
include { ASSESS_ASSEMBLY_COMPLEASM } from './assess_assembly_compleasm/main.nf'

workflow QC_WORKFLOW {
    take:
    ch_assembly

    main:
    ASSESS_ASSEMBLY_QUAST(ch_assembly)
    ASSESS_ASSEMBLY_COMPLEASM(ch_assembly)

    emit:
    quast = ASSESS_ASSEMBLY_QUAST.out
    compleasm = ASSESS_ASSEMBLY_COMPLEASM.out
}
