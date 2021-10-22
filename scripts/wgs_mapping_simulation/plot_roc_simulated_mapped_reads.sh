#!/bin/bash
# plot_roc_simulated_mapped_reads.sh: Preprocess the inputs for mapeval roc plots. Run the roc plot rendering script.

WORKDIR=${HOME}/run_sim_reads


# Where should temp files go?
mkdir -p "${WORKDIR}"
export TMPDIR="${WORKDIR}/tmp"
mkdir -p "${TMPDIR}"

cd ${WORKDIR}

# Replace all names of mappers with human-readable ones
function humanize_names() {
    sed -e 's/[a-zA-Z0-9_.]*bwamem[a-zA-Z0-9_.]*/BWA/' -e 's/[a-zA-Z0-9_.]*giraffe_parental_default[a-zA-Z0-9_.]*/GiraffeParental/' -e 's/[a-zA-Z0-9_.]*giraffe_snp1kg_default[a-zA-Z0-9_.]*/Giraffe1000GP/' -e 's/[a-zA-Z0-9_.]*giraffe_primary_default[a-zA-Z0-9_.]*/GiraffePrimary/'
}

declare -a REGION_LIST=( "high_conf_hg002_v4.2.1_regions_1M" "all_difficult_regions_hg002_v4.2.1_regions_1M" "alllowmapandsegdupregions_hg002_v4.2.1_regions_100M" "mhc_hg002_v4.2.1_regions_100M" "cmrg_hg002_v4.2.1_regions_100M" "high_conf_NO1000GP_hg002_v4.2.1_regions_1M" "all_difficult_regions_NO1000GP_hg002_v4.2.1_regions_1M" "alllowmapandsegdupregions_NO1000GP_hg002_v4.2.1_regions_100M" "mhc_hg002_NO1000GP_v4.2.1_regions_100M" "cmrg_hg002_NO1000GP_v4.2.1_regions_100M" )

for REGION in "${REGION_LIST[@]}" ; do
    cd ${WORK_DIR}
    GRAPH="hg002_sample_grch38"
    READS="HG002"
    PAIRING="paired"
    SPECIES="human"
    PE_OPTS="-- -pe"
    GBWT='sampled\.64'
    LINEAR_GRAPH="NOTAPPLICABLE"
    echo "Extracting toplot-${SPECIES}-headline_${GRAPH}-${READS}-${PAIRING}.tsv for region ${REGION}"
    cat bwamem_roc_stats.${REGION}.tsv | head -n1 > toplot-${SPECIES}-headline_${GRAPH}-${READS}-${PAIRING}.tsv

    # Grab linear BWA for this graph and these reads
    tail -q -n +2 bwamem_roc_stats.${REGION}.tsv | grep -P "(${GRAPH}(${GBWT})?${READS}|${LINEAR_GRAPH}(${GBWT})?${READS})" | grep ${PE_OPTS} | humanize_names >> toplot-${SPECIES}-headline_${GRAPH}-${READS}-${PAIRING}.tsv


    # Grab giraffe and map and graphaligner non-linear, and other linear mappers
    tail -q -n +2 giraffe_parental_roc_stats.${REGION}.tsv giraffe_snp1kg_roc_stats.${REGION}.tsv giraffe_primary_roc_stats.${REGION}.tsv | grep ${PE_OPTS} | grep -P "(${GRAPH}(${GBWT})?${READS}|${LINEAR_GRAPH}(${GBWT})?${READS})" | sed 's/null/0/g' | humanize_names >> toplot-${SPECIES}-headline_${GRAPH}-${READS}-${PAIRING}.tsv

    Rscript ${HOME}/vg_pedigree_manuscript_methods/wgs_mapping_simulation/plot-roc.R toplot-${SPECIES}-headline_${GRAPH}-${READS}-${PAIRING}.tsv roc-plot-${SPECIES}-headline_${GRAPH}-${READS}-${PAIRING}.${REGION}.png ${REGION}
done

