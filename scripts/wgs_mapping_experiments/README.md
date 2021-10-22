# Real Read Mapping Scripts

The scripts in this directory run the BWA-MEM, DRAGEN, VG Giraffe, and VG Pedigree mappers and workflows against real reads from the Genome-in-a-Bottle HG002 and HG005 samples.

## Running the scripts

First, run the bwa mem mapper on the HG002 and HG005 trio sample set.

```
./bwamem_map.sh
```

Next, run the Illumina DRAGEN platform mapper on the HG002 and HG005 trio sample set.

```
./dragen_map.sh
```

Following that, run the VG Giraffe alignment on the HG002 and HG005 trio sample sets against the 1000GP graph reference 

```
./giraffe_map.sh
```

Finally, run the VG Pedigree alignment workflow on the HG002 and HG005 trio sample sets.
Once, each, using the default deeptrio and deepvariant models.
And another, each using the trained deeptrio and deepvariant models.

```
./giraffe_pedigree_map.sh HG002 default
./giraffe_pedigree_map.sh HG005 default
./giraffe_pedigree_map.sh HG002 trained
./giraffe_pedigree_map.sh HG005 trained
```

## Output files

BWA-MEM alignment files will be located in `${HOME}/run_bwamem_mapping`.
DRAGEN alignment files will be located in `${HOME}/run_dragen_mapping`.
VG Giraffe alignments files will be located in `${HOME}/run_giraffe_mapping`.
VG Pedigree workflow alignment files will be located in `${HOME}/run_giraffe_pedigree_mapping`.

