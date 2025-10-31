process P_DUPS {
    tag "Purging haplotypic duplicates from merged $sample_id assembly with Purge_dups"
    publishDir "results/merged_assemblies", mode: 'copy'
    conda 'bioconda::purge_dups'

    input:
    tuple val(sample_id), val(assembler), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}_${assembler}_major_merged_purged.fa")

    stub:
    sample_id = sample_id
    assembler = assembler
    //run ONThill wrapper script here - 
    """
    touch "${sample_id}_${assembler}_major_merged_purged.fa"
    """
}
