process TRIM_READS_CHOPPER {
    tag "Trimming and filtering raw reads: $sample_id"
    cpus 10
    publishDir "results/reads/${trim_status}_reads", mode: 'copy'
    conda 'bioconda::nanofilt'

    input:
    tuple val(sample_id), val(trim_status), path(reads)

    output:
    tuple val(sample_id), val(trim_status), path("${sample_id}_${trim_status}.fq.gz")

    script:
    sample_id = sample_id
    reads = reads
    trim_status = 'trim'

    """
    gunzip -c $reads | \
    NanoFilt -q 10 -l 500 --headcrop 10 --tailcrop 10 | \
    gzip > ${sample_id}_${trim_status}.fq.gz
    """
}
