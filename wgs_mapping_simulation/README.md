# Simulated Read Mapping and Evaluation Scripts

These scripts runs sample graph construction, simulates reads from that graph, and runs mapping and evaluation of the simulated reads against the BWA-MEM aligner, the VG Giraffe graph aligner against the Primary and 1000GP graph and Parental graph.

## Running the scripts

First, run the construction workflow to generate the sample graph for read simulation.

```
./construct_sim_graph.v4.2.1_grch38.sh
```

Next, simulate reads from the high confident and difficult benchmark regions.

```
./sim_reads.sh
```

Following that, run the simulated mapping evaluation.
Mappers tested are BWA-MEM, VG Giraffe against the Primary graph, VG Giraffe against the 1000GP graph, and
VG Giraffe against the Parental graph as constructed by the VG Parental workflow using the default deepvariant models.

```
./run_mapevals.sh
```

Finally run the roc-plotting scripts to render mapping evaluation ROC curve plots.

```
./plot_roc_simulated_mapped_reads.sh
```

## Output files

All mapping output BAM and GAM files, simulated read FASTQ files, and evaluation output JSON and SVG files will be located in `${HOME}/run_sim_reads`.

