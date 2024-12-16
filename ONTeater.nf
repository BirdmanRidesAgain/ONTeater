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
    PILON as PILON_FLYE; PILON as PILON_ND; 
    QUAST as QUAST_FLYE; QUAST as QUAST_ND; 
    } from './modules/processes.nf'
include { QUICKMERGE; P_DUPS } from './modules/processes.nf' //merged assembly-related


workflow {
    rawreads_ch = Channel.fromList(get_name_file_pair(params.input_ONTreads))

    // Trim and visualize raw longread data
    NANOPLOT_RAW(rawreads_ch)
    trimreads_ch = NANOFILT(rawreads_ch) //trim
    trimreads_ch.view()
    NANOPLOT_TRIM(trimreads_ch)
    
    // Begin primary assembly
    trimreads_ch.view()
    pri_asm_flye_ch = FLYE(trimreads_ch)
        // Adds genome size estimate and core availability info to improve nextDenovo polishing
    nd_conf_ch = GET_NEXTDENOVO_PARAMS(pri_asm_flye_ch, params.nextdenovo_conf)
    pri_asm_nd_ch = NEXTDENOVO(trimreads_ch, nd_conf_ch)


    // Create channel of reads and assembled fastas for Racon polishing
    polish_flye_ch = pri_asm_flye_ch.join(trimreads_ch)
    //polish_nd_ch = pri_asm_nd_ch.join(trimreads_ch)
    
    // Polish and merge genomes
    racon_flye_ch = RACON_FLYE(polish_flye_ch)
    //racon_nd_ch = RACON_ND(polish_nd_ch)
    
    
    // Polish with Illumina data if Illumina data is non-null
        //Illumina reads MUST be paired and named like this: *{1,2}.fq.gz
        //Illumina sample IDs must also exactly match the longread sample names
    //shortreads_ch = Channel.fromFilePairs(params.input_shortreads)

    /*
    if (!shortreads_ch) { println("No shortreads found. Skipping Pilon polishing.")}
    else {
        println("Shortreads found: ")
        println("Polishing assemblies with Pilon")

        polish_flye_ch = PILON_FLYE(racon_flye_ch, shortreads_ch)
        polish_nd_ch = PILON_ND(racon_nd_ch, shortreads_ch)
    }
    */
    // Use quast to get n50 and number of large fragments for each assembly
    //quast_flye_ch = QUAST_FLYE(racon_flye_ch)
    //quast_nd_ch = QUAST_ND(racon_nd_ch)

    //merge and purge duplicate contigs
    //merged_ch = QUICKMERGE(quast_flye_ch, quast_nd_ch)
    //merged_purged_ch = P_DUPS(merged_ch) //replace p_dups call(s) with wrapper script.
    
    // QC AND VISUALIZE ASSEMBLED GENOME:
        /*depth_assess_ch = MOSDEPTH(merged_purged_ch)
        visual_output_ch = VISUALIZE(merged_purged_ch, depth_assess_ch)
         calls 'percent_of_genome_over_1mil.py'+'visualize_contig_lengths.R' to visualize dist. of contigs
         also calls 'flag_contig_depth.R' from mosdepth output - indicates probable mtDNA/bacterial contamination
         all called from ./modules/helper_scripts
         outputs .pdfs of relevant stats
         */
}