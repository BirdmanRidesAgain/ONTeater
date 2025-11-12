/*
 * Keiler Collier
 * ONTeater V1 - used to run genome assembly from concatenated, Kraken-filtered readfiles.
 *
 * Started: 21 Feb 2024
 * Last update: 12 Dec 2024
 * Pipeline input params supplied from nextflow.config
 * Invocation is nextflow run ONTeater.nf -stub -with-report -with-dag ONTeater.html
 */

log.info """\
    O N T E A T E R - N F   P I P E L I N E
    ===================================
    Project directory       : $projectDir
    Input ONT longreads     : ${params.input_ONTreads}
    Input PacBio longreads  : ${params.input_PBreads}
    Input shortreads        : ${params.input_shortreads}
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

//include { REMOVE_CONTAMINANTS } from './modules/remove_contaminants_kraken2/main.nf' //initial filtering for contaminants
include { VISUALIZE_READS_NANOPLOT as NANOPLOT_RAW; VISUALIZE_READS_NANOPLOT as NANOPLOT_TRIM } from './modules/visualize_reads_nanoplot/main.nf' //nanoplot-related
include { TRIM_READS_CHOPPER } from './modules/trim_reads_chopper/main.nf' //trimming
include { ASSEMBLE_FLYE } from './modules/assemble_flye/main.nf' //primary assemblers
include { GET_NEXTDENOVO_PARAMS } from './modules/get_nextdenovo_params/main.nf'
include { ASSEMBLE_NEXTDENOVO } from './modules/assemble_nextdenovo/main.nf'
include { //polishing-related
    POLISH_RACON as RACON_FLYE; POLISH_RACON as RACON_ND; 
    } from './modules/polish_racon/main.nf'


workflow {
    rawreads_ch = Channel.fromList(get_name_file_pair(params.input_ONTreads))

    // Trim and visualize raw longread data
    NANOPLOT_RAW(rawreads_ch)
    trimreads_ch = TRIM_READS_CHOPPER(rawreads_ch) //trim
    trimreads_ch.view()
    NANOPLOT_TRIM(trimreads_ch)
    
    // Begin primary assembly
    trimreads_ch.view()
    pri_asm_flye_ch = ASSEMBLE_FLYE(trimreads_ch)
        // Adds genome size estimate and core availability info to improve nextDenovo polishing
    nd_conf_ch = GET_NEXTDENOVO_PARAMS(pri_asm_flye_ch, params.nextdenovo_conf)
    pri_asm_nd_ch = ASSEMBLE_NEXTDENOVO(trimreads_ch, nd_conf_ch)


    // Create channel of reads and assembled fastas for Racon polishing
    polish_flye_ch = pri_asm_flye_ch.join(trimreads_ch)
    //polish_nd_ch = pri_asm_nd_ch.join(trimreads_ch)
    
    // Polish and merge genomes
    racon_flye_ch = RACON_FLYE(polish_flye_ch)
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
}