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

//include { REMOVE_CONTAMINANTS } from './modules/processes.nf' //initial filtering for contaminants
include { NANOPLOT as NANOPLOT_RAW; NANOPLOT as NANOPLOT_TRIM } from './modules/processes.nf' //nanoplot-related
include { NANOFILT; FLYE; GET_NEXTDENOVO_PARAMS; NEXTDENOVO} from './modules/processes.nf' //primary assemblers
include { //polishing-related
    RACON as RACON_FLYE; RACON as RACON_ND; 
    } from './modules/processes.nf'
include { QUAST_MERGE; QUICKMERGE } from './modules/processes.nf' //merged assembly-related
//include { P_DUPS } from './modules/processes.nf'

workflow {
    rawreads_ch = Channel.fromList(get_name_file_pair(params.input_ONTreads))

    // Trim and visualize raw longread data
    NANOPLOT_RAW(rawreads_ch)
    trimreads_ch = NANOFILT(rawreads_ch) //trim
    trimreads_ch.view()
    NANOPLOT_TRIM(trimreads_ch)
    trimreads_ch.view()

    // Begin primary assembly

    // Flye should run first - nondemanding assembler
    pri_asm_flye_ch = FLYE(trimreads_ch)
    polish_flye_ch = pri_asm_flye_ch.join(trimreads_ch)
    racon_flye_ch = RACON_FLYE(polish_flye_ch)

    // Then nextDenovo - requires parameters and must run without anything else
    nd_conf_ch = GET_NEXTDENOVO_PARAMS(params.nextdenovo_conf, params.genome_size) // Adds genome size estimate
    pri_asm_nd_ch = NEXTDENOVO(trimreads_ch, nd_conf_ch)
    polish_nd_ch = pri_asm_nd_ch.join(trimreads_ch)
    racon_nd_ch = RACON_ND(polish_nd_ch)
    
    //merge and purge duplicate contigs
    best_asm_ch = QUAST_MERGE(polish_flye_ch, polish_nd_ch)
    merged_ch = QUICKMERGE(polish_flye_ch, polish_nd_ch, best_asm_ch)
    
    // QC AND VISUALIZE ASSEMBLED GENOME:
        /*depth_assess_ch = MOSDEPTH(merged_purged_ch)
        visual_output_ch = VISUALIZE(merged_purged_ch, depth_assess_ch)
         calls 'percent_of_genome_over_1mil.py'+'visualize_contig_lengths.R' to visualize dist. of contigs
         also calls 'flag_contig_depth.R' from mosdepth output - indicates probable mtDNA/bacterial contamination
         all called from ./modules/helper_scripts
         outputs .pdfs of relevant stats
         */
}
