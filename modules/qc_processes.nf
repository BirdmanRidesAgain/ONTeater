// QC STEPS BELOW HERE
process QC_QUAST {
    tag "Assessing final $sample_id assembly contiguity with Quast"
    publishDir "results/${sample_id}/QC", mode: 'copy'
    conda 'bioconda::quast'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}_final_quast_report.txt")

    script:
    sample_id = sample_id

    """
    # run quast
    quast -ek $fasta --out ${sample_id}_final --no-html
    mv ${sample_id}_final/report.txt ${sample_id}_final_quast_report.txt
    """
}

process QC_COMPLEASM {
    tag "Assessing final $sample_id assembly completeness with compleasm"
    publishDir "results/${sample_id}/QC", mode: 'copy'
    conda 'bioconda::compleasm'

    input:
    tuple val(sample_id), path(fasta)
    val(BUSCO_lineage)

    output:
    tuple val(sample_id), path("${sample_id}_summary.txt")

    script:
    sample_id = sample_id
    lineage=BUSCO_lineage
    threads=40

    """
    # RUN COMPLEASM ON ASSEMBLY
    compleasm run -a $fasta -o ${sample_id}_compleasm -l $lineage -t $threads
    mv ${sample_id}_compleasm/summary.txt ${sample_id}_summary.txt
    """
}
