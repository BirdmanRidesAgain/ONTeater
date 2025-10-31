process FLYE {
    tag "Primary assembly of $sample_id with Flye"
    conda 'bioconda::flye'
    label 'parallel'

    input:
    tuple val(sample_id), val(trim_status), path(reads)

    output:
    tuple val(sample_id), val(assembler), path("${sample_id}_${assembler}.fa")

    script:
    sample_id = sample_id
    assembler = "Flye"
    num_iter = 3

    // params: threads are now set by the config file, not the process itself
    """
    # set variables:
    THREADS=\$((\$(nproc) / 2))

    echo "Primary assembly of $sample_id with $assembler"
    echo "\$THREADS allocated"

    # assemble with chosen assembler
    flye --nano-hq $reads --threads \$THREADS --iterations $num_iter -o ./
    mv assembly.fasta ${sample_id}_${assembler}.fa
    """
}
