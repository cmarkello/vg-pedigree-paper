#!/usr/bin/env bash
# This preprocesses the 1000 Genomes Project GRCh38 liftover VCFs to remove
# variants in large segmental duplications; not doing this produced a lot of
# miscalls in these regions realative to a graph without variants.

set -ex
set -o pipefail

function download() {
    if [ ! -e "${2}" ] ; then
        aws s3 cp --no-progress "${1}" "${2}"
    fi
}

function wget_download() {
    if [ ! -e "${2}" ] ; then
        wget_download "${1}" -O "${2}"
    fi
}

WORKDIR=${HOME}/run_graph_construction
WORKFLOW_INPUT_DIR=${WORKDIR}/inputs

# Where should temp files go?
mkdir -p "${WORKDIR}" "${WORKFLOW_INPUT_DIR}"
export TMPDIR="${WORKDIR}/tmp"
mkdir -p "${TMPDIR}"

cd ${WORKFLOW_INPUT_DIR}
rm -f GRCh38_segdups_gt10kb.bed.gz
wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/genome-stratifications/v2.0/GRCh38/SegmentalDuplications/GRCh38_segdups_gt10kb.bed.gz GRCh38_segdups_gt10kb.bed.gz
gzip -d GRCh38_segdups_gt10kb.bed.gz 
for CHROM in {1..22} X Y ; do
    (
        download s3://vg-data/1kg_GRCh38/variants/ALL.chr${CHROM}_GRCh38.genotypes.20170504.vcf.gz ALL.chr${CHROM}_GRCh38.genotypes.20170504.vcf.gz
        download s3://vg-data/1kg_GRCh38/variants/ALL.chr${CHROM}_GRCh38.genotypes.20170504.vcf.gz.tbi ALL.chr${CHROM}_GRCh38.genotypes.20170504.vcf.gz.tbi
        if [[ ! -e ALL.chr${CHROM}_GRCh38.genotypes.20170504.no_segdups_gt10kb.vcf.gz ]] ; then
            pbgzip -n 2 -dc ALL.chr${CHROM}_GRCh38.genotypes.20170504.vcf.gz | bcftools view -T ^GRCh38_segdups_gt10kb.bed -O v | pbgzip -n 2 -c > ALL.chr${CHROM}_GRCh38.genotypes.20170504.no_segdups_gt10kb.vcf.gz
        fi
        if [[ ! -e ALL.chr${CHROM}_GRCh38.genotypes.20170504.no_segdups_gt10kb.vcf.gz.tbi ]] ; then
            tabix -p vcf ALL.chr${CHROM}_GRCh38.genotypes.20170504.no_segdups_gt10kb.vcf.gz
        fi
    ) &
done

for CHROM in {1..22} X Y ; do
    wait
done

