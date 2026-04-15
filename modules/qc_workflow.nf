include { ASSESS_ASSEMBLY_QUAST } from './assess_assembly_quast/main.nf'
include { ASSESS_ASSEMBLY_COMPLEASM } from './assess_assembly_compleasm/main.nf'

workflow QC_WORKFLOW {
    take:
    ch_assembly

    main:
    ch_estimate_genome_bp_script = channel.value(file("${projectDir}/bin/estimate_genome_bp.py"))
    ASSESS_ASSEMBLY_QUAST(ch_assembly, ch_estimate_genome_bp_script)
    ASSESS_ASSEMBLY_COMPLEASM(ch_assembly)

    emit:
    quast = ASSESS_ASSEMBLY_QUAST.out
    compleasm = ASSESS_ASSEMBLY_COMPLEASM.out
}
