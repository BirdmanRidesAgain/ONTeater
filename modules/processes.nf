/*
 * KA Collier
 * ONTeater V1 - processes module
 * Started: Feb 21 2024
 * Last update: May 15 2024
 *
 * Helper script for ONTeater.nf, containing all processes.
 * 
 */

/*
process REMOVE_CONTAMINANTS {
    tag "Removing non-eukaryote reads: $sample_id"
}*/

process NANOPLOT {
    tag "Visualizing $trim_status reads: $sample_id"
    cpus 10
    publishDir "results/reads/read_stats/${trim_status}/", mode: 'copy'
    conda 'bioconda::nanoplot'

    input:
    tuple val(sample_id), val(trim_status), path(reads)

    output:
    path "${sample_id}_${trim_status}_NanoPlot"

    script:
    sample_id = sample_id
    trim_status = trim_status

    """
    echo "Running NanoPlot on $trim_status reads: $sample_id"
    echo "$task.cpus allocated"
    """
    if (trim_status == 'raw') {
        """
        NanoPlot -t $task.cpus --huge \
        --tsv_stats --drop_outliers --loglength --color limegreen \
        -p ${sample_id}_${trim_status}_ --title ${sample_id}_${trim_status} \
        --outdir ${sample_id}_${trim_status}_NanoPlot \
        --fastq $reads
        """
    } else {
        """
        NanoPlot -t $task.cpus --huge \
        --tsv_stats --drop_outliers --loglength --color royalblue \
        -p ${sample_id}_${trim_status}_ --title ${sample_id}_${trim_status} \
        --outdir ${sample_id}_${trim_status}_NanoPlot \
        --fastq $reads
        """ 
    }
}
process NANOFILT {
    tag "Trimming and filtering raw reads: $sample_id"
    cpus 10
    publishDir "results/reads/${trim_status}_reads", mode: 'copy'
    conda 'bioconda::nanofilt'

    input:
    tuple val(sample_id), val(trim_status), path(reads)

    output:
    tuple val(sample_id), val(trim_status), path("${sample_id}_${trim_status}.fq.gz")

    script:
    sample_id = sample_id
    reads = reads
    trim_status = 'trim'

    """
    gunzip -c $reads | \
    NanoFilt -q 10 -l 500 --headcrop 10 --tailcrop 10 | \
    gzip > ${sample_id}_${trim_status}.fq.gz
    """
}
process FLYE {
    tag "Primary assembly of $sample_id with Flye"
    conda 'bioconda::flye'
    label 'parallel'

    input:
    tuple val(sample_id), val(trim_status), path(reads)

    output:
    tuple val(sample_id), val(assembler), path("${sample_id}_${assembler}.fa")

    script:
    sample_id = sample_id
    assembler = "Flye"
    num_iter = 3

    // params: threads are now set by the config file, not the process itself
    """
    # set variables:
    THREADS=\$((\$(nproc) / 2))

    echo "Primary assembly of $sample_id with $assembler"
    echo "\$THREADS allocated"

    # assemble with chosen assembler
    flye --nano-hq $reads --threads \$THREADS --iterations $num_iter -o ./
    mv assembly.fasta ${sample_id}_${assembler}.fa
    """
}
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
process NEXTDENOVO {
    tag "Primary assembly of $sample_id with $assembler"
    conda 'bioconda::nextdenovo'
    label 'parallel'

    input:
    tuple val(sample_id), val(trim_status), path(reads)
    path "nextdenovo_conf"

    output:
    tuple val(sample_id), val(assembler), path("${sample_id}_${assembler}.fa")

    script: 
    sample_id = sample_id
    assembler = "nextDenovo"

    """
    # assemble with chosen assembler
    ls $reads > input_fofn
    nextDenovo $nextdenovo_conf
    ASSEMBLY=01_rundir/03.ctg_graph/nd.asm.fasta
    mv \$ASSEMBLY ${sample_id}_${assembler}.fa
    """
}

process RACON {
    tag "Polishing primary assembly $sample_id ($assembler) with Racon"
    publishDir "results/primary_assemblies/${assembler}", mode: 'copy'
    conda 'bioconda::racon bioconda::minimap2'
    label 'parallel'

    input:
    tuple val(sample_id), val(assembler), path(fasta), val(trim_status), path(reads)

    output:
    tuple val(sample_id), val(assembler), path("${sample_id}_${assembler}_racon.fa")

    script:
    Integer threads = 20
    sample_id = sample_id
    assembler = assembler

    """
    minimap2 -t $threads $fasta -ax map-ont $reads > ${assembler}.sam
    racon -t $threads -u $reads ${assembler}.sam $fasta > ${sample_id}_${assembler}_racon.fa
    """
}

/*
process QUICKMERGE {
    tag "Merging $sample_id assemblies with Quickmerge"
    publishDir "results/merged_assemblies", mode: 'copy'
    conda 'bioconda::quickmerge'

    input:
    tuple val(sample_id_1), val(assembler_1), val(n50_1), val(num_large_contigs_1), path(fasta_1) //this is flye
    tuple val(sample_id_2), val(assembler_2), val(n50_2), val(num_large_contigs_2), path(fasta_2) //this is nextDenovo

    output:
    tuple val(sample_id), val(assembler), path("${sample_id}_${assembler}_major_merged.fa")

    stub:
    sample_id = sample_id_1 //they're all the same sample_id, ffs
    if (n50_1 > n50_2) {
        println("$assembler_1 assembly more contiguous. Merging $fasta_2 into $fasta_1")
        assembler = assembler_1
        """
        touch "${sample_id}_${assembler}_major_merged.fa"
        """
    }
    else { //it's very unlikely that we'll have a tie.
        println("$assembler_2 assembly more contiguous. Merging $fasta_1 into $fasta_2")
        assembler = assembler_2
        """
        touch "${sample_id}_${assembler}_major_merged.fa"
        """ 
        }
    }
 
process P_DUPS {
    tag "Purging haplotypic duplicates from merged $sample_id assembly with Purge_dups"
    publishDir "results/merged_assemblies", mode: 'copy'
    conda 'bioconda::purge_dups'

    input:
    tuple val(sample_id), val(assembler), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}_${assembler}_major_merged_purged.fa")

    stub:
    sample_id = sample_id
    assembler = assembler
    //run ONThill wrapper script here - 
    """
    touch "${sample_id}_${assembler}_major_merged_purged.fa"
    """
}

/*
process MOSDEPTH {
    tag "Outputs depth-per-contig stats from input assembly"
}
*/

