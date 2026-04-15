process ASSEMBLE_NEXTDENOVO {
    tag "Primary assembly of $sample_id with $assembler"
    // nextDenovo currently breaks under Python 3.13 on this platform.
    // Pin Python to a known-compatible runtime for stable seed_cns execution.
    conda 'bioconda::nextdenovo=2.5.2 conda-forge::python=3.10'
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
