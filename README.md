# ONTeater
A NextFlow-enabled pipeline used to assemble genomes with ONT-only input. Illumina or Illumina-like (eg, BGI) shortreads are supported.

## Goals:
This is a pipeline meant to automate parallel assembly of genomes from ONT longread data. The modules and processes involved are based off my Genome_assembly repository.

## Installation:
The pipeline includes three dependencies: [NextFlow](https://www.nextflow.io/docs/latest/getstarted.html), and [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). You will need to install all three of these for the pipeline to run. [Docker](https://docs.docker.com/engine/install/) and [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) are not currently supported, but a Docker port is in the works.

## Workflow:
![DAG of ONTeater operation. Some functions are dry-run only at the moment.](https://github.com/BirdmanRidesAgain/ONTeater/blob/main/ONTeater_DAG.png?raw=true)
Figure 1. DAG graph of the pipeline execution.

Figure 2. Output from R visualization.

## Quickstart:
For a dry run with the included dummy data, the following code works:


nextflow run ONTeater.nf -stub -with-dag --input_rawreads toy_data/*ONT*fq.gz --input_shortreads toy_data/*{1,2}.fq.gz


By default, ONTeater will look for your input data files in the directory in which ONTeater.nf is launched - therefore, if you have multiple unique datasets to analyze, you can also just store them in the working directory. If you do this, be aware that the files must conform to a uniform naming convention (\*ONT\*fq.gz for long reads and *{1,2}.fq.gz for short reads). Otherwise, it's possible for the parser to mix up short and longread data, or fail to recognize a data file. *It's generally safer to input them as parameters.*

The pipeline was written with ONT data in mind, but will process PacBio data as well. PacBio data must be submitted directly as an argument to --input_rawreads, as the parser will not directly recognize them.


