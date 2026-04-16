# ONTeater

## Quickstart

A minimal example to assembly a human genome from ONT data can be run with the following command:

```sh
HUMAN_GENOME_SIZE=3.2g # also accepts g,m,k suffixes
BUSCO_LIN=primates
OUT_PREFIX=example_output
ONTeater --ONT_rds <input.fq.gz> --genome_size "$HUMAN_GENOME_SIZE" --BUSCO_lineage "$BUSCO_LIN" --out "$OUT_PREFIX"
```

Note that the user should supply a rough estimate of genome size and a valid [`BUSCO`](https://busco.ezlab.org) lineage.
If genome size is not roughly known, we recommend using the [GoaT](https://goat.genomehubs.org) database to estimate from a related species.  
[`NextFlow`](https://www.nextflow.io) and either [`conda`](https://anaconda.org/anaconda/conda) or [Docker](https://www.docker.com) are dependencies and assumed to be accessible inside your system.

Published Docker image (recommended):

```bash
docker pull ghcr.io/birdmanridesagain/onteater:1.0.0-amd64
```

## Table of Contents

- [Quickstart](#quickstart)
- [Introduction](#introduction)
- [High-level Description](#high-level-description)
- [Options](#options)

## Introduction

ONTeater is a [NextFlow](https://www.nextflow.io/docs/latest/index.html)-enabled pipeline behind a Python wrapper, intended to produce highly-contiguous genome assemblies with only ONT input.
Illumina shortreads and PacBio longreads are not currently supported, but may be in the future.

We initially developed the program for genome assembly in nonmodel organisms, with an emphasis on conservation.
ONTeater does not attempt to phase, but produces a single, collapsed haplotype.  

Thus, ONTeater performs well with relatively low-quality data, providing contiguous draft assemblies even for:

- Relatively short data
  - Input N50 of 5-10kbp
- Relatively small datasets
  - 10-20 Gbp

For extremely high-quality ONT datasets (ie, 100+ Gbp or input N50 >20kbp), we recommend [Hifiasm](https://github.com/chhylp123/hifiasm?tab=readme-ov-file#ontonly)'s ONT-only assembly mode, as this produces two phased assemblies.

## High-level Description

![ONTeater pipeline flowchart](./ONTeater_flowchart.png)

The `ONTeater` assembly pipeline takes ONT data as an input, optionally accepting PacBio and Illumina-like paired short reads as supplemental data sources.
It begins with data visualization and trimming with [`NanoPlot`](https://github.com/wdecoster/nanoplot) and [`Chopper`](https://github.com/wdecoster/chopper) (v.0.7.0; Wouter de Coster and Rademakers 2023), both from the `NanoPack` series of programs.
Low-quality reads (QUAL<10), reads under 500 base pairs in length, and ONT-specific adapters are all removed.

Trimmed data then are run through two different assembly algorithms – [`Flye`](https://github.com/mikolmogorov/Flye) (v.2.9.6; Kolmogorov et al. 2019, Lin et al. 2016), an older repeat-graph based assembler, and [`nextDenovo`]((https://github.com/nextomics/nextdenovo)) (v.2.5.2; Jiang et al. 2024), a string-graph based assembler.
These produce different graph walks (and hence, assemblies) throughout the same sets of reads, converging in areas of high complexity, but having variable performance around the graph ‘edges’.
These two ‘primary assemblies’ are individually polished with [`Racon`](https://github.com/isovic/racon) (v.1.5.0; Vaser et al. 2017).

The two assemblies are then merged using [`quickmerge`](https://github.com/mahulchak/quickmerge) (v.0.3; Chakraborty et al. 2016), using the least contiguous (as measured by N50) assembly to patch gaps in the more contiguous.

We then use [`purge_dups`](https://github.com/dfguan/purge_dups) (v.1.2.6; Guan et al. 2020) to attempt to collapse the assembly to a single haplotype.
Finally, a battery of QC tools, including but not limited to [`QUAST`](https://github.com/ablab/quast) (v.5.2.0; Gurevich et al. 2013), [`Compleasm`](https://github.com/huangnengCSU/compleasm) (v.0.2.6; Huang and Li 2023).
Results are then written to a final output directory for the end-user to consume.

## Options

All native options in `NextFlow` are usable in `ONTeater`.
Notably, `--trace` and `--report` are useful.

|Option|Default|Data type|Description|
|---|---|---|---|
|`--help`|NA|Flag|Set to print a help message and exit.|
|`--ONT_rds`|`null`|String|A path to a gzipped file of ONT reads used for assembly.|
|`--genome_size`|`null`|String|A value representing genome size. Bare numbers (eg, `3.2`) are interpreted as gigabasepairs (g); use `k`, `m`, or `g` suffixes for explicit units (eg, `250m`, `3.2g`).|
|`--BUSCO_lineage`|`null`|String|BUSCO lineage name for completeness assessment (see the [BUSCO lineage list](https://busco.ezlab.org/list_of_lineages.html)). Required for `run`, `trim`, `assemble`, and `qc`.|
|`--threads`|`null`|Integer|Thread count passed to Nextflow (`--threads`). If omitted, pipeline/default profile settings are used.|
|`--workflow`|`run`|String|Entry point for development/testing. Valid modes include: `run`, `trim`, `assemble`, `postprocess`, `qc`.|
|`--flye_asm`|`null`|String|A path to a `Flye` assembly. Used to bypass assembly and provide genomes directly to `postprocess` runs.|
|`--nd_asm`|`null`|String|A path to a `NextDenovo` assembly. Used to bypass assembly and provide genomes directly to `postprocess` runs.|
|`--final_asm`|`null`|String|Path to a final assembly FASTA used to run `qc` alone.|
|`--profile`|`null`|String|Nextflow profile(s) to use (eg, `standard`, `conda`, `docker`.|
|`--trace`|`null`|String|Enable Nextflow trace output and write to the provided file path (`-with-trace`).|
|`--report`|`null`|String|Enable Nextflow HTML report output and write to the provided file path (`-with-report`).|
|`--prefix`|`out`|String|Output prefix used under `results/<prefix>/...`.|
|`--container_image`|`ghcr.io/birdmanridesagain/onteater:1.0.0-amd64`|String|Docker image used with `-profile docker` (override for a custom image/tag). Most users will not need to edit this.|
|`--stub`|`false`|Flag|Run Nextflow native stub mode (`-stub-run`) to check dependencies on your system.|
|`--resume`|`false`|Flag|Resume a previous Nextflow run using cached work (`-resume`).|

## Execution profiles

### Conda (default)

Uses Nextflow’s built-in conda integration (`-profile standard` or `-profile conda`).
This works well in Linux environments, but has known issues in OSX, and has not been tested on Windows.

### Docker (recommended for reproducibility)

Prefer pulling a published tagged image, then running with the Docker profile:

```bash
docker pull ghcr.io/birdmanridesagain/onteater:1.0.0-amd64
nextflow run main.nf -profile docker --container_image ghcr.io/birdmanridesagain/onteater:1.0.0-amd64 --ONT_rds ... --genome_size ... --BUSCO_lineage ...
```

Wrapper equivalent:

```bash
ONTeater -pro docker --container_image ghcr.io/birdmanridesagain/onteater:1.0.0-amd64 --ONT_rds ... --genome_size ... --BUSCO_lineage ...
```

**Note:** Nextflow and your input data paths still run on the host; only **task processes** execute inside the container. Mounting and file permissions follow [Nextflow’s Docker documentation](https://www.nextflow.io/docs/latest/docker.html).

If you need to build locally instead (for development or custom changes):

```bash
docker buildx build --platform linux/amd64 -t ghcr.io/birdmanridesagain/onteater:1.0.0-amd64 -f docker/Dockerfile --load .
```

### Stub mode

Use Nextflow native stub mode for fast wiring/tests:

```bash
./ONTeater --workflow run --ONT_rds <reads.fq.gz> --genome_size 5m --BUSCO_lineage enterobacterales --stub
```

Equivalent native Nextflow call:

```bash
nextflow run main.nf --workflow run --ONT_rds <reads.fq.gz> --genome_size 5m --BUSCO_lineage enterobacterales -stub-run
```

## Smoke test

### Conda (default profile)

```bash
./scripts/smoke_test.sh
```

By default this uses **conda** (`standard` profile) and `--workflow run`, then validates outputs under `results/smoke_direct/` and `results/smoke_wrapper/` (including `reads/`, merged assemblies, and QC paths when applicable).

Useful overrides:

- `WORKFLOW_MODE=trim ./scripts/smoke_test.sh` for preprocess-only.
- `RUN_FULL=1 ./scripts/smoke_test.sh` to additionally launch a second full `workflow run` with prefix `smoke_full`.
- `WORKFLOW_MODE=postprocess FLYE_ASM=<flye.fa> ND_ASM=<nd.fa> ./scripts/smoke_test.sh`
- `WORKFLOW_MODE=qc FINAL_ASM=<final.fa> ./scripts/smoke_test.sh`

### Docker (`docker` profile)

Separate script so Docker and conda smoke paths do not interfere:

```bash
./scripts/smoke_test_docker.sh
```

This pulls/uses `ghcr.io/birdmanridesagain/onteater:1.0.0-amd64` by default, and builds from `docker/Dockerfile` only if the image is missing locally, then runs the same checks with `-profile docker` and distinct prefixes (`smoke_direct_docker`, `smoke_wrapper_docker`).

Optional:

- `CONTAINER_IMAGE=ghcr.io/birdmanridesagain/onteater:dev ./scripts/smoke_test_docker.sh`
- `SKIP_IMAGE_BUILD=1 ./scripts/smoke_test_docker.sh` if the image is already loaded locally.
- `STUB_RUN=1 ./scripts/smoke_test_docker.sh` to execute smoke with Nextflow native stubs.

### Requirements

The ONTeater pipeline is computationally heavy, but does not require GPU support to function.
The pipeline was primarily developed and tested on a remote OVH-BRUTE cluster with 755Gb of RAM and roughly 90 threads.
It is recommended that you run it on a cluster of similar or greater strength.

## Runtime

Overall runtime is strongly influenced by the target organism’s genome size and complexity, but scales roughly linearly with input data volume.
A representative mammalian assembly (*Stenella longirostris*, the spinner dolphin) at ~30x coverage was completed in roughly 66 hours.
More detailed statistics are in the work for this.
