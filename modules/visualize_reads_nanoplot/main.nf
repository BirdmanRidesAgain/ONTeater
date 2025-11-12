process VISUALIZE_READS_NANOPLOT {
    tag "Visualizing $trim_status reads: $sample_id"
    cpus 10
    publishDir "results/reads/read_stats/${trim_status}/", mode: 'copy'
    conda 'bioconda::nanoplot'

    input:
    tuple val(sample_id), val(trim_status), path(reads)

    output:
    path "${sample_id}_${trim_status}_NanoPlot"

    script:
    sample_id = sample_id
    trim_status = trim_status

    """
    echo "Running NanoPlot on $trim_status reads: $sample_id"
    echo "$task.cpus allocated"
    """
    if (trim_status == 'raw') {
        """
        NanoPlot -t $task.cpus --huge \
        --tsv_stats --drop_outliers --loglength --color limegreen \
        -p ${sample_id}_${trim_status}_ --title ${sample_id}_${trim_status} \
        --outdir ${sample_id}_${trim_status}_NanoPlot \
        --fastq $reads
        """
    } else {
        """
        NanoPlot -t $task.cpus --huge \
        --tsv_stats --drop_outliers --loglength --color royalblue \
        -p ${sample_id}_${trim_status}_ --title ${sample_id}_${trim_status} \
        --outdir ${sample_id}_${trim_status}_NanoPlot \
        --fastq $reads
        """ 
    }
}
