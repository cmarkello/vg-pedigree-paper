#!/usr/bin/env bash
# dragen_map.sh: map HG002, HG003, HG004, HG005, HG006, and HG007 150bp paired reads with Illumina's DRAGEN module to the HS38d1 reference

set -ex
set -o pipefail

function download() {
    if [ ! -e "${2}" ] ; then
        aws s3 cp --no-progress "${1}" "${2}"
    fi
}

function wget_download() {
    if [ ! -e "${2}" ] ; then
        wget "${1}" -O "${2}"
    fi
}

function copy() {
    if [ ! -e "${2}" ] ; then
        cp "${1}" "${2}"
    fi
}

SAMPLE_NAMES=("HG002" "HG003" "HG004" "HG005" "HG006" "HG007")
WORKDIR=${HOME}/run_dragen_mapping
REF_FASTA="GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna"

# Where should temp files go?
mkdir -p "${WORKDIR}"
export TMPDIR="${WORKDIR}/tmp"
mkdir -p "${TMPDIR}"

# Download data input data
cd $WORKDIR
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.gz "${WORKDIR}/${REF_FASTA}.gz"
gzip -d "${WORKDIR}/${REF_FASTA}.gz"

# Generate the DRAGEN reference index
mkdir -p ${WORKDIR}/dragen_index ${WORKDIR}/tmp && cd ${WORKDIR}
dragen --build-hash-table true \
  --output-directory ${WORKDIR}/dragen_index \
  --ht-reference ${WORKDIR}/${REF_FASTA}

for SAMPLE_NAME in "${SAMPLE_NAMES[@]}" ; do
    if [[ ${SAMPLE_NAME} == *"HG002"* ]]; then
        wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG002.novaseq.pcr-free.35x.R1.fastq.gz "${WORKDIR}/${SAMPLE_NAME}.R1.fastq.gz"
        wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG002.novaseq.pcr-free.35x.R2.fastq.gz "${WORKDIR}/${SAMPLE_NAME}.R2.fastq.gz"
    elif [[ ${SAMPLE_NAME} == *"HG003"* ]]; then
        wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG003.novaseq.pcr-free.35x.R1.fastq.gz "${WORKDIR}/${SAMPLE_NAME}.R1.fastq.gz"
        wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG003.novaseq.pcr-free.35x.R2.fastq.gz "${WORKDIR}/${SAMPLE_NAME}.R2.fastq.gz"
    elif [[ ${SAMPLE_NAME} == *"HG004"* ]]; then
        wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG004.novaseq.pcr-free.35x.R1.fastq.gz "${WORKDIR}/${SAMPLE_NAME}.R1.fastq.gz"
        wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG004.novaseq.pcr-free.35x.R2.fastq.gz "${WORKDIR}/${SAMPLE_NAME}.R2.fastq.gz"
    elif [[ ${SAMPLE_NAME} == *"HG005"* ]]; then
        wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG005.novaseq.pcr-free.30x.R1.fastq.gz "${WORKDIR}/${SAMPLE_NAME}.R1.fastq.gz"
        wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG005.novaseq.pcr-free.30x.R2.fastq.gz "${WORKDIR}/${SAMPLE_NAME}.R2.fastq.gz"
    elif [[ ${SAMPLE_NAME} == *"HG006"* ]]; then
        wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG006.novaseq.pcr-free.30x.R1.fastq.gz "${WORKDIR}/${SAMPLE_NAME}.R1.fastq.gz"
        wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG006.novaseq.pcr-free.30x.R2.fastq.gz "${WORKDIR}/${SAMPLE_NAME}.R2.fastq.gz"
    elif [[ ${SAMPLE_NAME} == *"HG007"* ]]; then
        wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG007.novaseq.pcr-free.30x.R1.fastq.gz "${WORKDIR}/${SAMPLE_NAME}.R1.fastq.gz"
        wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG007.novaseq.pcr-free.30x.R2.fastq.gz "${WORKDIR}/${SAMPLE_NAME}.R2.fastq.gz"
    fi
    
    # run bwa mem mapper
    cd $WORKDIR
    mkdir -p ${WORKDIR}/dragen_output_${SAMPLE_NAME}
    dragen -f \
      -r ${WORKDIR}/dragen_index \
      -1 "${SAMPLE_NAME}.R1.fastq.gz" \
      -2 "${SAMPLE_NAME}.R2.fastq.gz" \
      --RGID 1 \
      --RGSM ${SAMPLE_NAME} \
      --verbose --bin_memory=50000000000 --enable-map-align true --enable-variant-caller false \
      --pair-by-name=true \
      --enable-map-align-output=true \
      --intermediate-results-dir ${WORKDIR}/tmp \
      --output-directory ${WORKDIR}/dragen_output_${SAMPLE_NAME} \
      --output-file-prefix dragen_output_${SAMPLE_NAME} 2> dragen_output_${SAMPLE_NAME}.stderr
done



