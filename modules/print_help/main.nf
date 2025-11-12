def PRINT_HELP() {
    def message='''
    Basic Usage:
    ONTeater --ONT_rds path/to/variants.vcf --genome_size 1g --BUSCO_lineage <valid BUSCO lineage>
    For full parameter list, see README.md
    '''
    log.info message.stripIndent()
}