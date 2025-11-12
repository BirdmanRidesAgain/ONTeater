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
    // Get genome size from the Flye assembly if not defined in params
    Integer genome_size = params.genome_size
    if (params.genome_size == null) {
        genome_size = fasta.size()
    }
    """
    # Use sed to edit the appropriate line
    cat $nextdenovo_conf | sed "s/XXg/${genome_size}g/" > run.cfg
    """
}