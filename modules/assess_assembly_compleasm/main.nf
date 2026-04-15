process ASSESS_ASSEMBLY_COMPLEASM {
    tag "Assessing final $sample_id assembly completeness with compleasm"
    publishDir "results/QC", mode: 'copy'
    conda 'bioconda::compleasm bioconda::sepp'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}_final_compleasm_report.txt")

    script:
    sample_id = sample_id
    Integer threads = 40
    lineage = params.BUSCO_lineage

    if (lineage == null) { // use autolineage
        """
        compleasm run --autolineage -t $threads -a $fasta -o ${sample_id}
        mv ${sample_id}/summary.txt ${sample_id}_final_compleasm_report.txt
        """
    } else {
        """
        compleasm run -l $lineage -t $threads -a $fasta -o ${sample_id}
        mv ${sample_id}/summary.txt ${sample_id}_final_compleasm_report.txt
        """
    }

    stub:
    """
    printf "stub\\tASSESS_ASSEMBLY_COMPLEASM\\t%s\\n" "$sample_id" > ${sample_id}_final_compleasm_report.txt
    """
}
