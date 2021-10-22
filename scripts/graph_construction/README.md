# Graph Construction Scripts

These scripts run VG-based graph and index construction for the Primary and 1000GP-based graph references used in the subsequent mapping and variant calling experiments.

## Running the scripts

First, download the requisite vcf files and generate the segmental duplicate-free vcfs.

```
./preprocess_1000GPlons_GRCh38_graph.sh
```

Next, run the construction workflow for building the graph and indexes for the primary graph reference.

```
./make_primary_GRCh38_graph.sh
```

Finally, run the construction workflow for building the graph and indexes for the non-segmental-duplicate 1000GP-based graph reference.

```
./make_1000GPlons_GRCh38_graph.sh
```

## Output files

The Primary graph reference indexes will be located in `${HOME}/run_graph_construction/construct-grch38-primary-outstore`.
The 1000GP graph reference and indexes will be located in `${HOME}/run_graph_construction/construct-grch38-liftover-nosegdup-outstore`.

