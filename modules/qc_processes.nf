// QC STEPS BELOW HERE
process QC_QUAST {
    tag "Assessing final $sample_id assembly contiguity with Quast"
    publishDir "results/QC", mode: 'copy'
    conda 'bioconda::quast'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}_final_quast_report.txt")

    script:
    sample_id = sample_id
    boolean use_stub = params.qc_stub ?: false

    if (use_stub) {
        """
        printf "stub\\tQC_QUAST\\t%s\\n" "$sample_id" > ${sample_id}_final_quast_report.txt
        """
    } else {
        """
        # run quast
        quast -ek $fasta --out ${sample_id}_final --no-html
        mv ${sample_id}_final/report.txt ${sample_id}_final_quast_report.txt
        """
    }
}

process QC_COMPLEASM {
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
    boolean use_stub = params.qc_stub ?: false

    if (use_stub) {
        """
        printf "stub\\tQC_COMPLEASM\\t%s\\n" "$sample_id" > ${sample_id}_final_compleasm_report.txt
        """
    } else if (lineage == null) { //use autolineage
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
}
