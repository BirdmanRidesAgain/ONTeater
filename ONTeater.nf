/*
 * KA Collier
 * ONTeater V1 - used to QC fastas.
 *
 * Feb 16 2024
 * pipeline input params supplied from nextflow.config
 */

log.info """\
    O N T E A T E R - N F   P I P E L I N E
    ===================================
    Project directory   : $projectDir
    Input fasta         : ${params.input_fasta}
    Fasta QC directory  : ${params.fasta_QC_dir}
    """
    .stripIndent()

process QUAST {
    tag "Running Quast on $sample_id"
    publishDir params.fasta_QC_dir, mode: 'copy'
    //conda 'bioconda::quast'
    conda 'conda_envs/quast.yml'

    input:
    tuple val(sample_id), path(fasta)

    output:
    path "$sample_id"

    shell: // we're using 'shell' to allow bash (${}) vars AND nextflow vars(!{})
    '''
    quast -e -k !{fasta} -o !{sample_id}
    '''
}

def get_name_file_pair(paths) {
    // function to parse fasta paths and return labels for them as a tuple
    println "Parsing file names and locations: "
    path_list = files(paths) // creates list of filepaths
    int num_files = path_list.size

    def name_list = [] // creates empty list for the filenames
    path_list.each { name_list.add(it.getSimpleName()) } // adds each filename to a new list

    //now, we combine the lists with a for loop:
    def combined_list = []
    for (int i in 0..(num_files - 1)) {
        def fasta_tuple = new Tuple2(name_list[i], path_list[i])
        combined_list.add(fasta_tuple)
    }
    return combined_list
}

workflow {

    def fasta_list = get_name_file_pair(params.input_fasta)

    Channel 
        .fromList(fasta_list)
        .set { fasta_ch }
        //fasta_ch.subscribe onComplete: { println "List converted to channel"}
        fasta_ch.view()

    QUAST(fasta_ch)
}
