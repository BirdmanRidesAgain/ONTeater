process PURGE_HAPLOTYPIC_DUPLICATES {
    tag "Purging haplotypic duplicates from merged $sample_id assembly with Purge_dups"
    publishDir "results/merged_assemblies", mode: 'copy'

    input:
    tuple val(sample_id), val(assembler), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}_${assembler}_major_merged_purged.fa")

    script:
    //run ONThill wrapper script here - 
    """
    touch "${sample_id}_${assembler}_major_merged_purged.fa"
    """
    
    stub:
    //run ONThill wrapper script here - 
    """
    touch "${sample_id}_${assembler}_major_merged_purged.fa"
    """
}
