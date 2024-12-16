# ONTeater
A NextFlow-enabled pipeline used to assemble genomes with ONT-only input. Illumina or Illumina-like (eg, BGI) shortreads are supported.

## Notes specifically for Data Structures 2270:
Genome assembly was, frankly, a bad idea to implement for a final exam. This is because:
1. It requires some fun dependencies,
2. I can't include any non-trivial test sets, and 
3. Full runtime for most organisms takes 10-20 hours. 

However, DAGs are ubiquitous in bioinformatics, and this pipeline fit that criterion extremely well.

To get around problems 2-3, I have supplied both a stub function for every process in this pipeline, and 'toy' datasets, in the toy_data subdirectory. They have realistic filenames and will be picked up by the parser if input as arguments (or in the same directory as `nextflow.nf`), but are completely empty. The stub functions then make liberal use of *touch* to create further empty files, letting logic of the DAG (Figure 1) work without any of the computation, also defanging `NextFlow`'s tendency to crash HPC clusters in the process. Along the same lines, I commented out the checks for `Conda` and `Docker` environments, removing those as dependencies.

As discussed in QuickStart, this version of ONTeater can be run with:


`nextflow run ONTeater.nf -stub -with-dag --input_rawreads toy_data/*ONT*fq.gz --input_shortreads toy_data/*{1,2}.fq.gz`


Unfortunately, this still requires `NextFlow` (a weird derivative of Apache Groovy) itself as a dependency, and I assume that you probably don't want to install it. Therefore, I will include a quick explanation of what, exactly, is going on here.

Genome assembly is a problem that involved performing complicated, highly sequential manipulations on one or more input data files until you reach the stopping point of a finished genome. This is both inefficient and error-prone to do as a human, with many, many points at which mistakes can creep in. At the same time, there are enough decision points that a linear pipeline fails to accurately assemble. So, because there are no cycles in the assembly process (we don't go from finished products to raw data), it's modeled extremely well by a directed acyclic graph (DAG).

In `NextFlow`, you define 'processes' and 'channels', where the processes correspond to graph nodes and the channels are links between the nodes. The program is built around the idea of asynchronous parallel processing, so each 'process' in the DAG can be visited only one time. 

This creates substantial confusion in the DAG diagram itself - notice, for example, that several processes (NANOPLOT, RACON, PILON and QUAST) appear multiple times - for example, as 'NANOPLOT_RAW' and 'NANOPLOT_TRIM'. This makes no sense from a graphical perspective, but is a requirement in terms of `Nextflow`'s requirement for all process to be unique. Because, potentially, multiple datasets from different species can be run simultaneously (effectively multiplexing the 'channels'), ensuring that all nodes are visited only once ensures that data from different species do not collide.

This requirement can be hacked around, but it's not intended. An example of this contributes to an error seen in this DAG (regarding the PILON processes), which will be discussed later.

From a biological perspective, the interpretation of this graph is reasonably simple. Starting at the yellow box at the top, we read in our data, generating the first channels (Channel.fromList, Channel.fromFilePairs) - these represent our raw 'long-read' and 'short-read' data, respectively. The 'short-read' data is not necessarily present, and will be ignored if the channel is empty, resulting in a simpler graph.

We then clean and visualize our the long-read data (NANOPLOT and CHOPPER) - NANOPLOT produces quality-control graphs that the end-user wants to see, so its channel/output goes directly to the yellow box at the bottom. CHOPPER simply trims the data and passes the channel to the a pair of assembly algorithms (FLYE and NEXTDENOVO), splitting the DAG.

These assembly algorithms use further, much more complicated manipulations of de bruijin graphs in order to combine small (~10,000 base pairs) fragments of genomic information 'reads' to larger (~1,000,000 base pairs) ones 'contigs'. Because both implementations have somewhat different algorithms, we see differential performance between the two on data of different characteristics - FLYE works better on very repetitive data, and NEXTDENOVO performs better on shorter reads.

The channels continue to be split while errors (sequencing errors, misassemblies, etc...) are corrected using polishing algorithms. RACON uses mapping data from the assemblies themselves to correct, and PILON activates *only* if short reads were read taken in. Both processes then feed out to the same channel - however, because PILON's activation is governed by an if statement, it seems to be displaying oddly in the visual representation. (Either that, or I've made an error, which is also very possible.)

We then assess which draft assembly is higher quality and use that to guide a merge, patching breaks in the assembly (QUAST, QUICKMERGE). There's a final, post-processing step to remove suspected duplicate contigs.

Further steps are in various levels of completeness, but all would fit into the same general 'expand onto the graph' framework.


## Goals:
This is a pipeline meant to automate parallel assembly of genomes from ONT longread data. The modules and processes involved are based off my Genome_assembly repository.

## Installation:
The pipeline includes three dependencies: [NextFlow](https://www.nextflow.io/docs/latest/getstarted.html), and [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). You will need to install all three of these for the pipeline to run. [Docker](https://docs.docker.com/engine/install/) and [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) are not currently supported, but a `Docker` port is in the works.

## Workflow:
![DAG of ONTeater operation. Some functions are dry-run only at the moment.](https://github.com/BirdmanRidesAgain/ONTeater/blob/main/ONTeater_DAG.png?raw=true)
Figure 1. DAG graph of the pipeline execution.

![Sample output from 'visualize_contig_lengths.R' module.](https://github.com/BirdmanRidesAgain/ONTeater/blob/main/modules/visualize_descending_contigs_sample_output.png?raw=true)
Figure 2. Sample output from `visualize_contig_lengths.R` module. Each bar represents a contig, with red bars denoting contigs over 1,000,000 bp in length. Green bars (not shown in this assembly) represent smaller ones, with the horizontal red line denoting the 1-million bp mark.

## Quickstart:
For a dry run with the included dummy data, the following code works:


`nextflow run ONTeater.nf -stub -with-dag --input_rawreads toy_data/*ONT*fq.gz --input_shortreads toy_data/*{1,2}.fq.gz`


By default, ONTeater will look for your input data files in the directory in which ONTeater.nf is launched - therefore, if you have multiple unique datasets to analyze, you can also just store them in the working directory. If you do this, be aware that the files must conform to a uniform naming convention (\*ONT\*fq.gz for long reads and *{1,2}.fq.gz for short reads). Otherwise, it's possible for the parser to mix up short and longread data, or fail to recognize a data file. *It's generally safer to input them as parameters.*

The pipeline was written with ONT data in mind, but will process PacBio data as well. PacBio data must be submitted directly as an argument to --input_rawreads, as the parser will not directly recognize them.


