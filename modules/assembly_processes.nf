/*
 * KA Collier
 * ONTeater V1 - processes module
 * Started: Feb 21 2024
 * Last update: Feb 4 2025
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
    //tool = "nanoplot" not actually needed here but included for consistency
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
    tool = "nanofilt"
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
    THREADS=\$(nproc)

    echo "Primary assembly of $sample_id with $assembler"
    echo "\$THREADS allocated"

    # assemble with chosen assembler
    flye --nano-hq $reads --threads \$THREADS --iterations $num_iter -o ./
    mv assembly.fasta ${sample_id}_${assembler}.fa
    """
}
process GET_NEXTDENOVO_PARAMS { // this can also be implemented as a helper function for NEXTDENOVO itself.
    tag "Fetching params for NextDenovo assembly of $sample_id"
    publishDir 'results'

    input:
    //tuple val(sample_id), val(assembler), path(fasta)
    path(nextdenovo_conf)
    val(genome_size)

    output:
    path("nd_run_goodparams.cfg")


    script:
        // Get genome size from the Flye assembly if not defined in params
    //not functional yet
    /*
    if (params.genome_size == null) {
        genome_size = fasta.size()
    }
    */
    Integer genome_size = genome_size
    """
    # Use sed to edit the appropriate line
    cat $nextdenovo_conf | sed "s/XXg/${genome_size}g/" > nd_run_goodparams.cfg
    """
}
process NEXTDENOVO {
    tag "Primary assembly of $sample_id with $assembler"
    conda 'bioconda::nextdenovo'
    label 'parallel'
    maxForks 1
    cpus Runtime.runtime.availableProcessors()
    
    input:
    tuple val(sample_id), val(trim_status), path(reads)
    path (nextdenovo_conf)

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
    Integer threads = 80
    sample_id = sample_id
    tool = "racon"
    assembler = assembler

    """
    minimap2 -t $threads $fasta -ax map-ont $reads > ${assembler}.sam
    racon -t $threads -u $reads ${assembler}.sam $fasta > ${sample_id}_${assembler}_racon.fa
    """
}

process QUAST_MERGE {
    tag "Determining which $sample_id assembly is most contiguous (Quast)"
    publishDir "results/primary_assemblies", mode: 'copy'
    conda 'bioconda::quast'

    input:
    tuple val(sample_id_1), val(assembler_1), path(fasta_1) //This is Flye
    tuple val(sample_id_2), val(assembler_2), path(fasta_2) //This is nextDenovo

    output:
    tuple val(sample_id), path("best_N50_asm.txt")

    script:
    sample_id = sample_id_1 //FIXME - there needs to be a check here to ensure that ID1 and ID2 are identical

    """
    # run quast
    quast -ek $fasta_1 --out ${sample_id_1}_${assembler_1} --no-html
    quast -ek $fasta_2 --out ${sample_id_2}_${assembler_2} --no-html
    ASM1_N50=\$(cat ${sample_id_1}_${assembler_1}/report.txt | grep "N50" | awk '{print \$NF}')
    ASM2_N50=\$(cat ${sample_id_2}_${assembler_2}/report.txt | grep "N50" | awk '{print \$NF}')

    # Compare N50s and output best assembler to textfile
        # If both are equal, select ASM2 (conventionally nextDenovo)
    if [ \$ASM1_N50 -gt \$ASM2_N50 ]
    then
        echo $assembler_1 > best_N50_asm.txt
    else
        echo $assembler_2 > best_N50_asm.txt
    fi        
    """
}


process QUICKMERGE {
    tag "Merging $sample_id assemblies with Quickmerge"
    publishDir "results/merged_assemblies", mode: 'copy'
    conda 'bioconda::quickmerge'

    input:
    tuple val(sample_id_1), val(assembler_1), path(fasta_1) //this is flye
    tuple val(sample_id_2), val(assembler_2), path(fasta_2) //this is nextDenovo
    tuple val(sample_id), path(best_N50_asm) //this is a textfile bc restrictions with channels

    output:
    tuple val(sample_id), path("merged_${sample_id_1}_${assembler_1}_major.fasta")
    //tuple val(sample_id), val(assembler), path("${sample_id}_${assembler}_major_merged.fa")

    script:
    if (sample_id_1 == sample_id_2) //checks to see that the assemblies we're merging are the same animal.
        sample_id = sample_id_1
        """
        ASM_MAJOR=\$(cat $best_N50_asm)
        if [ \$ASM_MAJOR == $assembler_1 ]
        then
            merge_wrapper.py -pre ${sample_id_1}_${assembler_1}_major $fasta_1 $fasta_2
        else
            merge_wrapper.py -pre ${sample_id_2}_${assembler_2}_major $fasta_2 $fasta_1
        fi
        """
        //find a way to raise an error if this fucks up       
}

process P_DUPS {
    tag "Purging haplotypic duplicates from merged $sample_id assembly with purge_dups"
    publishDir "results/merged_assemblies", mode: 'copy'
    conda 'bioconda::purge_dups'

    input:
    tuple val(sample_id_fa), path(fasta)
    tuple val(sample_id_rds), val(trim_status), path(reads) // we do not need to have the sample ID here; the second one isn't used. Just here b/c easier to have it.

    output:
    tuple val(sample_id), path("${sample_id}_merged.purged.fa")

    script:
    sample_id = sample_id_rds
    Integer threads = 80

    //run ONTeater wrapper script here - 
    """
    INPUT_PAF=${sample_id}.paf.gz
    INPUT_SPLIT=${sample_id}.split
    SPLIT_PAF=${sample_id}.split.self.paf.gz

    minimap2 -t $threads -x map-ont $fasta $reads | pigz > \$INPUT_PAF
    pbcstat \$INPUT_PAF #produces PB.stat and PB.base.cov
    calcuts PB.stat > cutoffs 2> calcuts.log
    split_fa $fasta > \$INPUT_SPLIT
    minimap2 -t $threads -x asm5 -DP \$INPUT_SPLIT \$INPUT_SPLIT | pigz -c > \$SPLIT_PAF

    purge_dups -2 -T cutoffs -c PB.base.cov \$SPLIT_PAF > dups.bed 2> purge_dups.log
    get_seqs dups.bed -e $fasta -p ${sample_id}_merged
    """
}

