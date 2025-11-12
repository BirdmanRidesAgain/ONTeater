#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
/*
 * Keiler Collier
 * ONTeater V1 - used to run genome assembly from concatenated, Kraken-filtered readfiles.
 *
 * Started: 21 Feb 2024
 * Last update: 12 Nov 2025
 */

log.info """\
    O N T E A T E R - N F   P I P E L I N E
    ===================================
    Project directory       : $projectDir
    Input ONT longreads     : ${params.ONT_rds}
    Input PacBio longreads  : ${params.PB_rds}
    Input shortreads        : ${params.short_rds}
    Genome size             : ${params.genome_size}
    BUSCO lineage           : ${params.BUSCO_lineage}
    """
    .stripIndent()

def get_name_file_pair(paths) {
    // function to parse fasta paths and return labels for them as a tuple
    path_list = files(paths) // creates list of filepaths
    int num_files = path_list.size

    def name_list = [] // creates empty list for the filenames
    path_list.each { name_list.add(it.getSimpleName()) } // adds each filename to a new list

    //now, we combine the lists with a for loop:
    def combined_list = []
    for (int i in 0..(num_files - 1)) {
        def file_tuple = new Tuple3(name_list[i], 'raw', path_list[i])
        combined_list.add(file_tuple)
    }
    return combined_list
}

def get_sample_id(String path) {
    //Helper function to provide a channel of sample names. Avoids tuple issues when only running parts of pipeline.
    def infile = file(path)
    return infile.getSimpleName()
}
//include { REMOVE_CONTAMINANTS } from './modules/remove_contaminants_kraken2/main.nf' //initial filtering for contaminants
// Import modules
include { PRINT_HELP } from './modules/print_help'
include { VISUALIZE_READS_NANOPLOT as NANOPLOT_RAW; VISUALIZE_READS_NANOPLOT as NANOPLOT_TRIM } from './modules/visualize_reads_nanoplot/main.nf' //nanoplot-related
include { TRIM_READS_CHOPPER } from './modules/trim_reads_chopper/main.nf' //trimming
include { ASSEMBLE_FLYE } from './modules/assemble_flye/main.nf' //primary assemblers
include { GET_NEXTDENOVO_PARAMS } from './modules/get_nextdenovo_params/main.nf'
include { ASSEMBLE_NEXTDENOVO } from './modules/assemble_nextdenovo/main.nf'
include { //polishing-related
    POLISH_RACON as RACON_FLYE; POLISH_RACON as RACON_ND; 
    } from './modules/polish_racon/main.nf'

def print_help() {
    //Prints out ONTeater instructions.
    log.info"""
    Basic Usage:
    nextflow run ONTeater.nf [--ONT_raw OR --PB_raw]
    
    Options:
        --help          Flag. Show this help message and exit
        --ONT_raw       String; default null. A gzipped file of ONT reads used for assembly. Can be used with --PB_raw.
        --PB_raw        String; default null. A gzipped file of PacBio reads used for assembly. Can be used with --ONT_raw.
        --genome_size   Float; default 1.0. Genome size in GB
        --BUSCO_lineage String; default null (ie, compleasm auto-lineage). See https://busco.ezlab.org/list_of_lineages.html for acceptable list
        --workflow      String; default 'run'. Determines start point of workflow; valid options are: 'run', 'trim', 'assemble', 'merge', 'pdups', 'qc'
        --flye_asm      String; default null. Used in 'merge', 'pdups' and 'qc' workflows.
        --nd_asm        String; default null. Used in 'merge', 'pdups' and 'qc' workflows.
        --trace         Flag. Nextflow-native; produces a trace of the run.
        --report        String. Nextflow-native; produces a runtime report with the supplied name.

    Notes:
        The only mandatory option is either --ONT_raw or --PB_raw. The assembler also accepts both at the same time.
        Any errors (of which there are probably many) should be reported to https://github.com/BirdmanRidesAgain.
    
        The 'workflow' parameter is likely to be particularly buggy- this was a feature implemented for testing purposes.
    """.stripIndent()
}

//include { REMOVE_CONTAMINANTS } from './modules/processes.nf' //initial filtering for contaminants
include { NANOPLOT as NANOPLOT_RAW; NANOPLOT as NANOPLOT_TRIM; NANOFILT;  //nanoplot-related
    FLYE; GET_NEXTDENOVO_PARAMS; NEXTDENOVO;  //primary assemblers
    RACON as RACON_FLYE; RACON as RACON_ND;  //polishing-related
    QUAST_MERGE; QUICKMERGE; P_DUPS;  //merged assembly-related
} from './modules/assembly_processes.nf'
include { QC_QUAST; QC_COMPLEASM } from './modules/qc_processes.nf' //merged assembly-related

workflow {
    main:

    if (params.help) {
        PRINT_HELP()
        exit 10
    }

    rawreads_ch = Channel.fromList(get_name_file_pair(params.ONT_rds))

    // Trim and visualize raw longread data
    NANOPLOT_RAW(rawreads_ch)
    trimreads_ch = TRIM_READS_CHOPPER(rawreads_ch) //trim
    trimreads_ch.view()
    NANOPLOT_TRIM(trimreads_ch)
    
    // Begin primary assembly
    trimreads_ch.view()
    //pri_asm_flye_ch = ASSEMBLE_FLYE(trimreads_ch)
        // Adds genome size estimate and core availability info to improve nextDenovo polishing
    //nd_conf_ch = GET_NEXTDENOVO_PARAMS(pri_asm_flye_ch, params.nextdenovo_conf)
    //pri_asm_nd_ch = ASSEMBLE_NEXTDENOVO(trimreads_ch, nd_conf_ch)


    // Create channel of reads and assembled fastas for Racon polishing
    //polish_flye_ch = pri_asm_flye_ch.join(trimreads_ch)
    //polish_nd_ch = pri_asm_nd_ch.join(trimreads_ch)
    
    // Polish and merge genomes
    //racon_flye_ch = RACON_FLYE(polish_flye_ch)
    //racon_nd_ch = RACON_ND(polish_nd_ch)
    
    //merge and purge duplicate contigs
    //merged_ch = MERGE_ASSEMBLIES_QUICKMERGE(quast_flye_ch, quast_nd_ch)
    //merged_purged_ch = PURGE_HAPLOTYPIC_DUPLICATES(merged_ch) //replace p_dups call(s) with wrapper script.
    
    // QC AND VISUALIZE ASSEMBLED GENOME:
        /*depth_assess_ch = MOSDEPTH(merged_purged_ch)
        visual_output_ch = VISUALIZE(merged_purged_ch, depth_assess_ch)
         calls 'percent_of_genome_over_1mil.py'+'visualize_contig_lengths.R' to visualize dist. of contigs
         also calls 'flag_contig_depth.R' from mosdepth output - indicates probable mtDNA/bacterial contamination
         all called from ./modules/helper_scripts
         outputs .pdfs of relevant stats
         */
    
//    publish:
//    final = Channel.("Hello World")
}

//output {
//    final {}
//}
