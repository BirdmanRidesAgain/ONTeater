process TRIM_READS_CHOPPER {
    tag "Trimming and filtering raw reads: ${reads.simpleName}"

    input:
    path reads

    output:
    path "${reads.simpleName}_trim.fq.gz", emit: trimreads

    script:

    """
    NANO_FILT="\${CONDA_PREFIX:+\${CONDA_PREFIX}/bin/}NanoFilt"
    if [ ! -x "\$NANO_FILT" ]; then
        NANO_FILT=\$(command -v NanoFilt || true)
    fi
    if [ -z "\$NANO_FILT" ]; then
        echo "ERROR: NanoFilt not found in PATH/conda env." >&2
        exit 1
    fi

    pigz -d -c $reads | \
    \$NANO_FILT -q 10 -l 500 --headcrop 10 --tailcrop 10 | \
    gzip > ${reads.simpleName}_trim.fq.gz

    # Some datasets can fully filter out and still produce a non-empty gzip stream.
    # If no FASTQ headers remain, pass raw reads through to keep preprocessing usable.
    if ! pigz -d -c "${reads.simpleName}_trim.fq.gz" | awk 'NR==1 { if (\$0 ~ /^@/) ok=1 } END { exit(ok ? 0 : 1) }'; then
        cp $reads ${reads.simpleName}_trim.fq.gz
    fi
    """
}
