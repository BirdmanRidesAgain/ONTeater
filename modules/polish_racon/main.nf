process RACON {
    tag "Polishing primary assembly $sample_id ($assembler) with Racon"
    publishDir "results/primary_assemblies/${assembler}", mode: 'copy'
    conda 'bioconda::racon bioconda::minimap2'
    label 'parallel'

    input:
    tuple val(sample_id), val(assembler), path(fasta), val(trim_status), path(reads)

    output:
    tuple val(sample_id), val(assembler), path("${sample_id}_${assembler}_racon.fa")

    script:
    Integer threads = 20
    sample_id = sample_id
    assembler = assembler

    """
    minimap2 -t $threads $fasta -ax map-ont $reads > ${assembler}.sam
    racon -t $threads -u $reads ${assembler}.sam $fasta > ${sample_id}_${assembler}_racon.fa
    """
}
