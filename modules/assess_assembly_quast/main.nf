process ASSESS_ASSEMBLY_QUAST {
    tag "Assessing final $sample_id assembly contiguity with Quast"
    publishDir "results/QC", mode: 'copy'
    conda 'bioconda::quast'

    input:
    tuple val(sample_id), path(fasta)
    path estimate_genome_bp_script

    output:
    tuple val(sample_id), path("${sample_id}_final_quast_report.txt"), path("${sample_id}_final_quast_report.tsv"), path("${sample_id}_final_quast_report.pdf")

    script:
    sample_id = sample_id

    """
    set -euo pipefail

    estimated_bp=\$(python3 $estimate_genome_bp_script --user-size "${params.genome_size ?: ''}" --flye-asm "${params.flye_asm ?: ''}" --assembly "$fasta")

    quast_extra=""
    if [[ "\$estimated_bp" -gt 100000000 ]]; then
      quast_extra="--large -k"
    fi

    quast.py $fasta \$quast_extra -o ${sample_id}_final
    mv ${sample_id}_final/report.txt ${sample_id}_final_quast_report.txt
    mv ${sample_id}_final/report.tsv ${sample_id}_final_quast_report.tsv
    mv ${sample_id}_final/report.pdf ${sample_id}_final_quast_report.pdf
    """

    stub:
    """
    printf "stub\\tASSESS_ASSEMBLY_QUAST\\t%s\\n" "$sample_id" > ${sample_id}_final_quast_report.txt
    printf "stub\\tASSESS_ASSEMBLY_QUAST\\t%s\\n" "$sample_id" > ${sample_id}_final_quast_report.tsv
    touch ${sample_id}_final_quast_report.pdf
    """
}
