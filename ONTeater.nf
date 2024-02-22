/*
 * Keiler Collier
 * ONTeater V1 - used to run genome assembly from concatenated, Kraken-filtered readfiles.
 *
 * Feb 16 2024
 * Pipeline input params supplied from nextflow.config
 */

log.info """\
    O N T E A T E R - N F   P I P E L I N E
    ===================================
    Project directory   : $projectDir
    Input longreads     : ${params.input_rawreads}
    Input shortreads    : ${params.input_shortreads}
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

include { NANOPLOT as NANOPLOT_RAW; NANOPLOT as NANOPLOT_TRIM } from './modules/processes.nf'
include { CHOPPER; FLYE; GET_NEXTDENOVO_PARAMS; NEXTDENOVO} from './modules/processes.nf'
include {RACON as RACON_FLYE; RACON as RACON_ND; ASSESS as ASSESS_FLYE; ASSESS as ASSESS_ND } from './modules/processes.nf'
workflow {
    rawreads_ch = Channel.fromList(get_name_file_pair(params.input_rawreads))

    // Trim and visualize raw longread data
    NANOPLOT_RAW(rawreads_ch)
    trimreads_ch = CHOPPER(rawreads_ch) //trim
    NANOPLOT_TRIM(trimreads_ch)
    
    // Begin primary assembly
    pri_asm_flye_ch = FLYE(trimreads_ch)
        // Adds genome size estimate and core availability info to improve nextDenovo polishing
    nd_conf_ch = GET_NEXTDENOVO_PARAMS(pri_asm_flye_ch, params.nextdenovo_conf)
    pri_asm_nd_ch = NEXTDENOVO(trimreads_ch, nd_conf_ch)


    // Create channel of reads and assembled fastas for Racon polishing
    polish_flye_ch = pri_asm_flye_ch.join(trimreads_ch)
    polish_nd_ch = pri_asm_nd_ch.join(trimreads_ch)
    
    // Polish and merge genomes
    racon_flye_ch = RACON_FLYE(polish_flye_ch)
    racon_nd_ch = RACON_ND(polish_nd_ch)
    // Polish with Illumina data if Illumina data is non-null
    
    racon_flye_ch.view()
    assess_flye_ch = ASSESS_FLYE(racon_flye_ch)
    /*
    if (params.input_shortreads != null) {
        println "Shortreads found at params.shortreads"
        println "Polishing assemblies with Pilon"
        pilon_flye_ch = PILON_FLYE(racon_flye_ch, params.shortreads)
        pilon_nd_ch = PILON_ND(racon_nd_ch, params.shortreads)

        assess_flye_ch = ASSESS_FLYE(pilon_flye_ch)
        assess_nd_ch = ASSESS_ND(pilon_nd_ch)
    } else {
        println "no shortreads found"
        assess_flye_ch = ASSESS_FLYE(racon_flye_ch)
        assess_nd_ch = ASSESS_ND(racon_nd_ch)
    }
    */
    //assess_nd_ch.view()
    //merged_ch = QUICKMERGE(racon_flye_ch, racon_nd_ch)


    // Code for running QC on a fasta.
    //QUAST(fasta_ch)
    
}
