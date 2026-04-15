process GET_NEXTDENOVO_PARAMS { // this can also be implemented as a helper function for NEXTDENOVO itself.
    // We need to edit this command so that it also changes the nextDenovo resource use.
        // Otherwise, it's going to crash the Omen and other resource-constrained computers
    tag "Fetching params for NextDenovo assembly of $sample_id"
    publishDir 'results'

    input:
    tuple val(sample_id), val(assembler), path(fasta)
    path nextdenovo_conf

    output:
    path "run.cfg"

    script:
    // Prefer explicit genome size (e.g. 5m, 1.2g). If absent, estimate from Flye assembly size.
    String genome_size = params.genome_size ? params.genome_size.toString() : String.format('%.3fg', fasta.size() / 1e9)
    """
    # Replace the template token with concrete genome size.
    sed "s/__GENOME_SIZE__/${genome_size}/" $nextdenovo_conf > run.cfg
    """
}