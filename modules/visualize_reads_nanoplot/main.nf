process VISUALIZE_READS_NANOPLOT {
    tag "Visualizing $trim_status reads: ${reads.simpleName}"

    input:
    path reads
    val trim_status

    output:
    path "${reads.simpleName}_${trim_status}_NanoPlot", emit: viz

    script:

    if (trim_status == 'raw') {
        color = 'limegreen'
    } else { color = 'royalblue' }
        """
        NANO_PLOT="\${CONDA_PREFIX:+\${CONDA_PREFIX}/bin/}NanoPlot"
        if [ ! -x "\$NANO_PLOT" ]; then
            NANO_PLOT=\$(command -v NanoPlot || true)
        fi
        if [ -z "\$NANO_PLOT" ]; then
            echo "ERROR: NanoPlot not found in PATH/conda env." >&2
            exit 1
        fi
        \$NANO_PLOT -t $task.cpus --huge \
        --tsv_stats --drop_outliers --loglength --color ${color} \
        -p ${reads.simpleName}_${trim_status} --title ${reads.simpleName}_${trim_status} \
        --outdir ${reads.simpleName}_${trim_status}_NanoPlot \
        --fastq $reads
        """
    }
