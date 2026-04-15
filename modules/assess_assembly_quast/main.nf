process ASSESS_ASSEMBLY_QUAST {
    tag "Assessing final $sample_id assembly contiguity with Quast"
    publishDir "results/QC", mode: 'copy'
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

    stub:
    """
    printf "stub\\tASSESS_ASSEMBLY_QUAST\\t%s\\n" "$sample_id" > ${sample_id}_final_quast_report.txt
    """
}
