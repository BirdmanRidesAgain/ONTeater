/*
process REMOVE_CONTAMINANTS {
    tag "Removing non-eukaryote reads: $sample_id"
}*/

/*

/*
process MOSDEPTH {
    tag "Outputs depth-per-contig stats from input assembly"
}
*/


/*
process P_DUPS {
    tag "Purging haplotypic duplicates from merged $sample_id assembly with Purge_dups"
    publishDir "results/merged_assemblies", mode: 'copy'
    conda 'bioconda::purge_dups'

    input:
    tuple val(sample_id), path(fasta)

    output:
    val(sample_id)

    script:
    Integer threads = 80
    //run ONTeater wrapper script here - 
    """
    minimap2 -t $threads -x map-ont
    """
}
*/