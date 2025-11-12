process ASSEMBLE_NEXTDENOVO {
    tag "Primary assembly of $sample_id with $assembler"
    conda 'bioconda::nextdenovo'
    label 'parallel'

    input:
    tuple val(sample_id), val(trim_status), path(reads)
    path nextdenovo_conf

    output:
    tuple val(sample_id), val(assembler), path("${sample_id}_${assembler}.fa")

    script: 
    sample_id = sample_id
    assembler = "nextDenovo"

    """
    # assemble with chosen assembler
    ls $reads > input_fofn
    nextDenovo $nextdenovo_conf
    ASSEMBLY=01_rundir/03.ctg_graph/nd.asm.fasta
    mv \$ASSEMBLY ${sample_id}_${assembler}.fa
    """
}
