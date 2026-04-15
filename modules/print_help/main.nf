def PRINT_HELP() {
    def message='''
    Basic Usage:
    ONTeater --ONT_rds tests/data/barcode42_ONT_salmonella.fq.gz --genome_size 5m --BUSCO_lineage enterobacterales --prefix test_run

    Native Nextflow usage:
    nextflow run main.nf --ONT_rds tests/data/barcode42_ONT_salmonella.fq.gz --genome_size 5m --BUSCO_lineage enterobacterales --prefix test_run

    Notes:
    --ONT_rds is the ONT input flag.
    --workflow modes:
      run        preprocess + primary assembly + postprocess + qc
      trim       preprocess only
      assemble   primary assembly only
      postprocess requires --flye_asm and --nd_asm
      qc          requires --final_asm
    For full parameter list, see README.md
    '''
    log.info message.stripIndent()
}