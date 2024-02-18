# ONTeater
A NextFlow-enabled pipeline used to assemble genomes with ONT-only input.

## Goals:
This is a pipeline meant to automate parallel assembly of genomes from ONT longread data. The modules and processes involved are based off my Genome_assembly repository.

## Installation:
The pipeline includes three dependencies: [NextFlow](https://www.nextflow.io/docs/latest/getstarted.html), [Docker](https://docs.docker.com/engine/install/), and [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). You will need to install all three of these for the pipeline to run. [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) is not currently supported, but ports are welcome.

## Workflow:
Figure 1. DAG graph of the pipeline execution.

## Quickstart:
nextflow ONTeater.nf -resume
