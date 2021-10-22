# VG Pedigree Mapper, Simulation, and Variant Calling Evaluation Scripts

This repository contains scripts used to reproduce our work with the vg pedigree workflow in [toil-vg](https://github.com/vgteam/toil-vg) to conduct the mapping, read simulation, and variant calling experiments.

## Workflow Overview

The scripts should be run in the following order:

* [Data preparation and graph construction workflow scripts](scripts/graph_construction)
* [Read simulation and simulated mapping evaluation scripts](scripts/wgs_mapping_simulation)
* [Mapping scripts, for assessing accuracy between VG Pedigree alignments and other competing mappers](scripts/wgs_mapping_experiments)
* [Variant calling scripts and evaluation, for assessing variant calling accuracy from alignments produced between VG Pedigree and other competing mappers](scripts/wgs_calling_experiments)

## Output files

All output files by default will locate under directories placed within the path pointed to by the `${HOME}` environment variable.
Graph construction outputs will be located within `${HOME}/run_graph_construction`.
Simulated reads and mapping evaluations will be located within `${HOME}/run_sim_reads`.
Mapping output from BWAMEM alignments will be located within `${HOME}/run_bwamem_mapping`.
Mapping output from DRAGEN alignments will be located within `${HOME}/run_dragen_mapping`.
Mapping output from VG Giraffe alignments will be located within `${HOME}/run_giraffe_mapping`.
Mapping output from VG Pedigree alignments will be located within `${HOME}/run_giraffe_pedigree_mapping`.
Calling output and evaluation data from DeepTrio and Dragen variant callers will be located within `${HOME}/run_genotyping`.

## Replication Considerations

Though the top level workflows and scripts are automated, they have a number of system-level requirements. Access to a server that houses the Illumina Dragen Bio-IT Platform is required to run the DRAGEN-dependent scripts in `scripts/wgs_mapping_experiments` and `scripts/wgs_calling_experiments`.
The `toil-vg` program used in `scripts/graph_construction`, `scripts/wgs_mapping_simulation`, and `scripts/wgs_mapping_experiments` requires a significant amount of memory and scratch disk space to run in single-machine mode that these scripts are configured to run on by default. A server with at least a 32 core CPU, 100GB of memory and at least 5TB of local disk space is recommended to run in single-machine mode.
All software in these scripts rely on software containers in order to run, local installation of Docker is required to run them.



