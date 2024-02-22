/*
 * KA Collier
 * ONTeater V1 - processes module
 * Feb 21 2024
 *
 * Helper script for ONTeater.nf, containing all processes.
 * 
 */

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

    if (trim_status == 'raw') {
        """
        echo "Running NanoPlot on reads: $sample_id"
        echo "$task.cpus allocated"

        NanoPlot -t $task.cpus --huge \
        --tsv_stats --drop_outliers --loglength --color limegreen \
        -p ${sample_id}_${trim_status}_ --title ${sample_id}_${trim_status} \
        --outdir ${sample_id}_${trim_status}_NanoPlot \
        --fastq $reads
        """
    } else {
        """
        echo "Running NanoPlot on reads: $sample_id"
        echo "$task.cpus allocated"

        NanoPlot -t $task.cpus --huge \
        --tsv_stats --drop_outliers --loglength --color royalblue \
        -p ${sample_id}_${trim_status}_ --title ${sample_id}_${trim_status} \
        --outdir ${sample_id}_${trim_status}_NanoPlot \
        --fastq $reads
        """ 
    }
}
process CHOPPER {
    tag "Trimming and filtering raw reads: $sample_id"
    cpus 10
    publishDir "results/reads/${trim_status}_reads", mode: 'copy'
    conda 'bioconda::chopper'

    input:
    tuple val(sample_id), val(trim_status), path(reads)

    output:
    tuple val(sample_id), val(trim_status), path("${sample_id}_${trim_status}.fq.gz")

    script:
    sample_id = sample_id
    trim_status = 'trim'

    """
    gunzip -c $reads | \
    chopper --threads $task.cpus \
    -q 10 -l 500 --headcrop 10 --tailcrop 10 | \
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
    //threads = <figure out how to define threads in groovy>

    """
    # set variables:
    THREADS=\$((\$(nproc) / 2)) # 50% of your available processors. We go hard.

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
    cat $nextdenovo_conf | sed "s/genome_size = ...g/genome_size = $genome_size/" > run.cfg
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
process ASSESS {
    // This runs Quast on an assembly and returns the assembly+metadata, plus num_ctgs and N50
    tag "Calculating $sample_id ($assembler) N50 and number of contigs with Quast"
    publishDir "results/primary_assemblies/${assembler}", mode: 'copy'
    conda 'bioconda::quast'

    input:
    tuple val(sample_id), val(assembler), path(fasta)

    output:
    tuple val(sample_id), val(assembler), val(n50), val(num_large_contigs), path(fasta)

    exec:
    sample_id = sample_id
    assembler = assembler
    fasta = fasta
    n50 = '0'
    num_large_contigs = '0'
    
    println "Calculating $sample_id ($assembler) N50 and number of contigs with Quast"
    //"""
   // quast -e -k ${fasta} -o ${sample_id}_quast
   // echo "${sample_id}_quast/report.txt"
   // ${num_large_contigs}=`cat ${sample_id}_quast/report.txt | grep "# contigs (>= 50000 bp)" | awk '{ print \$NF }'`
   // ${n50}=`cat ${sample_id}_quast/report.txt | grep "N50" | awk '{ print \$NF }'`
   // """
    println "Assembly $sample_id ($assembler) N50: $n50"
    println "Assembly $sample_id ($assembler) number large (>50,000bp) contigs: $num_large_contigs"
}
/*
process QUICKMERGE {
    tag "Merging assemblies with Quickmerge"
    publishDir "results/merged_assemblies", mode: 'copy'
    conda 'bioconda::quickmerge'

    input:
    tuple val(sample_id1), val(assembler1), Integer(n50_1), Integer(num_large_contigs1), path(fasta1)
    tuple val(sample_id2), val(assembler2), Integer(n50_2), Integer(num_large_contigs2), path(fasta2)

    output:
    tuple val(sample_id3), path(fasta3)

    script:
    if (num_large_contigs1 >= num_large_contigs2 && n50_1 >= n50_2) {
        println "$assembler1 assembly more contiguous. Merging $assembler2 into $assembler1"
    } 
    if else (num_large_contigs2 >= num_large_contigs1 && n50_1 >= n50_2)
}

process QUAST {
    tag "Running Quast on $sample_id"
    publishDir params.fasta_QC_dir, mode: 'copy'
    conda 'bioconda::quast'

    input:
    tuple val(sample_id), path(fasta)

    output:
    path "quast/${sample_id}"

    script: 
    """
    quast -e -k ${fasta} -o quast/${sample_id}
    """

}
*/
