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
    Project directory           : $projectDir
    Selected workflow           : ${params.workflow}
    Input ONT longreads         : ${params.ONT_raw}
    Input PacBio longreads      : ${params.PB_raw}
    Input shortreads            : ${params.shortreads_raw}
    Input Flye assembly         : ${params.flye_asm}
    Input nextDenovo assembly   : ${params.nd_asm}
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

def get_sample_id(paths) {
    //Helper function to provide a channel of sample names. Avoids tuple issues when only running parts of pipeline.
    path_list = files(paths)
    int num_files = path_list.size

    def name_list = []
    path_list.each { name_list.add(it.getSimpleName()) } // adds each filename to a new list
    return name_list
}

//include { REMOVE_CONTAMINANTS } from './modules/processes.nf' //initial filtering for contaminants
include { NANOPLOT as NANOPLOT_RAW; NANOPLOT as NANOPLOT_TRIM } from './modules/processes.nf' //nanoplot-related
include { NANOFILT; FLYE; GET_NEXTDENOVO_PARAMS; NEXTDENOVO} from './modules/processes.nf' //primary assemblers
include { RACON as RACON_FLYE; RACON as RACON_ND; } from './modules/processes.nf' //polishing-related
include { QUAST_MERGE; QUICKMERGE; P_DUPS } from './modules/processes.nf' //merged assembly-related
//include { P_DUPS } from './modules/processes.nf'

workflow {
// INTRODUCTORY BEHAVIOR
    //CHECK INPUT PARAMETERS:
    //almost everything requires reads
    if (params.workflow != 'qc') {
        if (params.ONT_raw == null & params.PB_raw == null) { 
            error "Long reads required for 'run', 'trim', 'assemble', 'merge' and 'pdups' modes. Set --ONT_raw or --PB_raw." 
            }
    }

    // FIXME - this needs to become an 'input' channel, which can then be edited to process different read combos
    rawreads_ch = Channel.from(get_name_file_pair(params.ONT_raw))


    //This is really only used when you're running 'merge', 'pdups' or 'qc' as solo items
    def sample_id 
    sample_id = get_sample_id(params.ONT_raw)


    // Full complement of processes- desired behavior when workflow='run'
    DO_TRIM=true; DO_ASSEMBLE=true; DO_MERGE = true; DO_P_DUPS = true; DO_QC = true

    // This block checks what processes and inputs we need based on the workflow
    if (params.workflow != "run") {
        if (params.workflow == 'trim') {
            println "Identical workflow to 'run"
        } else if (params.workflow == "assemble") {
            DO_TRIM=false;
        } else if (params.workflow == "merge") {
            if (params.flye_asm == null | params.nd_asm == null) {
                error "Assemblies must be supplied for 'merge' workflow. Set --flye_asm and --nd_asm"
            }
            DO_TRIM=false; DO_ASSEMBLE=false;
        } else if (params.workflow == "pdups") {
            if (params.flye_asm == null | params.nd_asm == null) {
                error "Assemblies must be supplied for 'pdups' workflow. Set --flye_asm and --nd_asm."
            }
            DO_TRIM=false; DO_ASSEMBLE=false; DO_MERGE = false;
        } else if (params.workflow == 'qc') {
            if (params.flye_asm == null | params.nd_asm == null) {
                error "Assemblies must be supplied for 'qc' workflow. Set --flye_asm and --nd_asm."
            }
            DO_TRIM=false; DO_ASSEMBLE=false; DO_MERGE = false; DO_P_DUPS = false;
        } else { 
            println "Invalid workflow." 
            println "Valid options are: 'run', 'trim', 'assemble', 'merge', 'pdups', 'qc'."
        }
    }

// BEGIN BIOINFORMATIC PROCESSING
    if (DO_TRIM) {
        // Trim and visualize raw longread data
        NANOPLOT_RAW(rawreads_ch)
        trimreads_ch = NANOFILT(rawreads_ch)
        //trimreads_ch.view()
        NANOPLOT_TRIM(trimreads_ch)
    } else {
        //recreate the assembly channels from parameter input
        trimreads_ch = rawreads_ch
    }

    if (DO_ASSEMBLE) {
        // Begin primary assembly
        // Flye should run first - nondemanding assembler
        pri_asm_flye_ch = FLYE(trimreads_ch)
        polish_flye_ch = pri_asm_flye_ch.join(trimreads_ch)
        racon_flye_ch = RACON_FLYE(polish_flye_ch)
        // Then nextDenovo - requires parameters and must run without anything else
        nd_conf_ch = GET_NEXTDENOVO_PARAMS(params.nextdenovo_conf, params.genome_size) // Adds genome size estimate //FIXME - this calls for sampleID and sampleID not given
        pri_asm_nd_ch = NEXTDENOVO(trimreads_ch, nd_conf_ch)
        polish_nd_ch = pri_asm_nd_ch.join(trimreads_ch)
        racon_nd_ch = RACON_ND(polish_nd_ch)
    }
    else {
        flye_tup = new Tuple3 ("$sample_id", 'Flye', files(params.flye_asm))
        polish_flye_ch = Channel.of(flye_tup)
        nd_tup = new Tuple3 ("$sample_id", 'nextDenovo', files(params.nd_asm))
        polish_nd_ch = Channel.of(nd_tup)
        } 

    if (DO_MERGE) {
        //merge and purge duplicate contigs  
        best_asm_ch = QUAST_MERGE(polish_flye_ch, polish_nd_ch)
        merge_ch = QUICKMERGE(polish_flye_ch, polish_nd_ch, best_asm_ch)
    } else {
        //FIXME - make it so that you don't have to set 'flye_asm'; it can be any asm.
        merge_tup = new Tuple2 ("$sample_id", files(params.flye_asm))
        merge_ch = Channel.of(merge_tup)
    }

    if (DO_P_DUPS) {
        merge_purge_ch = P_DUPS(merge_ch, trimreads_ch)
    }

    /*
    if (DO_QC) {
        // QC AND VISUALIZE ASSEMBLED GENOME:
        /*depth_assess_ch = MOSDEPTH(merged_purged_ch)
        visual_output_ch = VISUALIZE(merged_purged_ch, depth_assess_ch)
        calls 'percent_of_genome_over_1mil.py'+'visualize_contig_lengths.R' to visualize dist. of contigs
        also calls 'flag_contig_depth.R' from mosdepth output - indicates probable mtDNA/bacterial contamination
        all called from ./modules/helper_scripts
        outputs .pdfs of relevant stats
    }
    */
}
