include { MERGE_ASSEMBLIES_QUICKMERGE } from './merge_assemblies_quickmerge/main.nf'
include { PURGE_HAPLOTYPIC_DUPLICATES } from './purge_haplotypic_duplicates/main.nf'

process GET_ASSEMBLY_CONTIGUITY {
    tag "Collecting contiguity stats for $sample_id ($assembler)"

    input:
    tuple val(sample_id), val(assembler), path(fasta)
    path contiguity_script

    output:
    tuple val(sample_id), val(assembler), path("n50.txt"), path("num_large_contigs.txt"), path("${sample_id}_${assembler}.fa")

    script:
    """
    cp $fasta ${sample_id}_${assembler}.fa
    python3 $contiguity_script --fasta ${sample_id}_${assembler}.fa --n50-out n50.txt --num-large-out num_large_contigs.txt
    """

    stub:
    """
    cp $fasta ${sample_id}_${assembler}.fa
    echo "0" > n50.txt
    echo "0" > num_large_contigs.txt
    """
}

workflow POSTPROCESS_ASSEMBLY {
    take:
    ch_flye_polished
    ch_nextdenovo_polished

    main:
    ch_primary_polished = ch_flye_polished.mix(ch_nextdenovo_polished)
    ch_contiguity_script = channel.value(file("${projectDir}/bin/calc_contiguity.py"))
    GET_ASSEMBLY_CONTIGUITY(ch_primary_polished, ch_contiguity_script)
    ch_contiguity_stats = GET_ASSEMBLY_CONTIGUITY.out.map { sample_id, assembler, n50_file, num_large_contigs_file, fasta ->
        def n50 = n50_file.text.trim().isLong() ? (n50_file.text.trim() as Long) : 0L
        def num_large_contigs = num_large_contigs_file.text.trim().isInteger() ? (num_large_contigs_file.text.trim() as Integer) : 0
        tuple(sample_id, assembler, n50, num_large_contigs, fasta)
    }
    ch_flye_for_merge = ch_contiguity_stats.filter { _sample_id, assembler, _n50, _num_large_contigs, _fasta -> assembler == 'Flye' }
    ch_nd_for_merge = ch_contiguity_stats.filter { _sample_id, assembler, _n50, _num_large_contigs, _fasta -> assembler == 'nextDenovo' }

    MERGE_ASSEMBLIES_QUICKMERGE(ch_flye_for_merge, ch_nd_for_merge)
    ch_merged = MERGE_ASSEMBLIES_QUICKMERGE.out

    PURGE_HAPLOTYPIC_DUPLICATES(ch_merged)

    emit:
    merged_assembly = ch_merged
    purged_assembly = PURGE_HAPLOTYPIC_DUPLICATES.out.purged_assembly
}
