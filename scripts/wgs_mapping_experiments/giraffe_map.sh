#!/usr/bin/env bash
# giraffe_map.sh: map HG002, HG003, 150bp and 250bp paired reads with VG GIRAFFE, VG GIRAFFE FAST, and VG GIRAFFE PRIMARY to the HS38d1-based graph references

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
WORKDIR=${HOME}/run_giraffe_mapping
REF_FASTA="GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna"

# Where should temp files go?
mkdir -p "${WORKDIR}"
export TMPDIR="${WORKDIR}/tmp"
mkdir -p "${TMPDIR}"

# Download data input data
cd $WORKDIR
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/liftover_snp1kg_grch38_nosegdup.xg -O ${WORKDIR}/snp1kg_decoys.xg
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/liftover_snp1kg_grch38_nosegdup.gbwt -O ${WORKDIR}/snp1kg_decoys.gbwt
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/liftover_snp1kg_grch38_nosegdup.gg -O ${WORKDIR}/snp1kg_decoys.gg
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/liftover_snp1kg_grch38_nosegdup.min -O ${WORKDIR}/snp1kg_decoys.min
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/liftover_snp1kg_grch38_nosegdup.dist -O ${WORKDIR}/snp1kg_decoys.dist
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.dict -O ${WORKDIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.dict

for SAMPLE_NAME in "${SAMPLE_NAMES[@]}" ; do
    SEQ_DICT="GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.dict"
    XG_IDX="-x snp1kg_decoys.xg"
    GBWT_IDX="-H snp1kg_decoys.gbwt"
    MIN_IDX="-m snp1kg_decoys.min"
    GG_IDX="-g snp1kg_decoys.gg"
    DIST_IDX="-d snp1kg_decoys.dist"
    FAST_PARAM=""
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
    
    # run giraffe mapper
    cd $WORKDIR
    docker run \
    -e SEQ_DICT=${SEQ_DICT} \
    -e XG_IDX=${XG_IDX} \
    -e GBWT_IDX=${GBWT_IDX} \
    -e MIN_IDX=${MIN_IDX} \
    -e GG_IDX=${GG_IDX} \
    -e DIST_IDX=${DIST_IDX} \
    -e READ1="${SAMPLE_NAME}.R1.fastq.gz" \
    -e READ2="${SAMPLE_NAME}.R2.fastq.gz" \
    -v ${PWD}:${HOME} -w ${HOME} quay.io/vgteam/vg:v1.31.0 \
    vg giraffe \
    -t 16 \
    ${XG_IDX} \
    ${GBWT_IDX} \
    ${MIN_IDX} \
    ${GG_IDX} \
    ${DIST_IDX} \
    -f ${READ1} \
    -f ${READ2} \
    ${FAST_PARAM} \
    --output-format BAM \
    --ref-paths ${SEQ_DICT} > giraffe_${SAMPLE_NAME}.bam
done

