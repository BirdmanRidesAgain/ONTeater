process ASSESS_ASSEMBLY_COMPLEASM {
    tag "Assessing final $sample_id assembly completeness with compleasm"
    conda 'bioconda::compleasm bioconda::sepp'

    input:
    tuple val(sample_id), path(fasta)
    path run_compleasm_script

    output:
    tuple val(sample_id), path("${sample_id}_final_compleasm_report.txt")

    script:
    sample_id = sample_id
    Integer threads = 40
    lineage = params.BUSCO_lineage ?: ''

    """
    python3 $run_compleasm_script \
      --assembly $fasta \
      --sample-id ${sample_id} \
      --threads $threads \
      --lineage "$lineage" \
      --report-out ${sample_id}_final_compleasm_report.txt
    """

    stub:
    """
    printf "stub\\tASSESS_ASSEMBLY_COMPLEASM\\t%s\\n" "$sample_id" > ${sample_id}_final_compleasm_report.txt
    """
}
