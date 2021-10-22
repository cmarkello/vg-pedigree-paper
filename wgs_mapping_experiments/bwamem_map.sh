#!/usr/bin/env bash
# bwamem_map.sh: map HG002, HG003, 150bp and 250bp paired reads with BWAMEM to the HS38d1 reference

set -ex
set -o pipefail

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
WORKDIR=${HOME}/run_bwamem_mapping
REF_FASTA="GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna"

# Where should temp files go?
mkdir -p "${WORKDIR}"
export TMPDIR="${WORKDIR}/tmp"
mkdir -p "${TMPDIR}"

# Download data input data
cd $WORKDIR
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.gz "${WORKDIR}/${REF_FASTA}.gz"
gzip -d "${WORKDIR}/${REF_FASTA}.gz"

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
    docker run \
    -e REF_FASTA=${REF_FASTA} \
    -e READ1="${SAMPLE_NAME}.R1.fastq.gz" \
    -e READ2="${SAMPLE_NAME}.R2.fastq.gz" \
    -v ${PWD}:${HOME} -w ${HOME} biocontainers/bwa:v0.7.17_cv1 \
    bwa mem \
    -t 16 \
    ${REF_FASTA} \
    ${READ1} \
    ${READ2} \
    > bwamem_${SAMPLE_NAME}.sam
    
done

