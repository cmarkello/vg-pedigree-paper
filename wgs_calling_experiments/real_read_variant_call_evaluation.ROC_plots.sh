#!/usr/bin/env bash
# real_read_variant_call_evaluation.ROC_plots.sh: Run rtg rocplot on rtg vcfeval output to generate ROC plots.

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

function run_default_rtgrocplot() {
    if [ ! -e "${1}.svg" ] ; then
        OUT_NAME=${1}
        TITLE="${2}"
        ZOOM_RANGE_X_START=${3}
        ZOOM_RANGE_X_END=${4}
        ZOOM_RANGE_Y_START=${5}
        ZOOM_RANGE_Y_END=${6}
        DRAGEN_RTG_ROC=${7}
        BWAMEM_RTG_ROC=${8}
        VG_1000GP_DEFAULT_RTG_ROC=${9}
        VG_PARENTAL_DEFAULT_RTG_ROC=${10}
        docker run \
        -e OUT_NAME=${OUT_NAME} \
        -e TITLE="${TITLE}" \
        -e ZOOM_RANGE_X_START=${ZOOM_RANGE_X_START} \
        -e ZOOM_RANGE_X_END=${ZOOM_RANGE_X_END} \
        -e ZOOM_RANGE_Y_START=${ZOOM_RANGE_Y_START} \
        -e ZOOM_RANGE_Y_END=${ZOOM_RANGE_Y_END} \
        -e DRAGEN_RTG_ROC=${DRAGEN_RTG_ROC} \
        -e BWAMEM_RTG_ROC=${BWAMEM_RTG_ROC} \
        -e VG_1000GP_DEFAULT_RTG_ROC=${VG_1000GP_DEFAULT_RTG_ROC} \
        -e VG_PARENTAL_DEFAULT_RTG_ROC=${VG_PARENTAL_DEFAULT_RTG_ROC} \
        -v ${PWD}:${HOME} -w ${HOME} realtimegenomics/rtg-tools:3.8.4 \
            rtg rocplot \
            --title "${TITLE}" \
            --plain \
            --zoom=${ZOOM_RANGE_X_START},${ZOOM_RANGE_Y_START},${ZOOM_RANGE_X_END},${ZOOM_RANGE_Y_END} \
            --curve=${BWAMEM_RTG_ROC}=bwamem \
            --curve=${DRAGEN_RTG_ROC}=dragen \
            --curve=${VG_1000GP_DEFAULT_RTG_ROC}=vg_1000gp \
            --curve=${VG_PARENTAL_DEFAULT_RTG_ROC}=vg_parent \
            --svg=${OUT_NAME}.svg
    fi
}

function run_trained_chr20_rtgrocplot() {
    if [ ! -e "${1}.svg" ] ; then
        OUT_NAME=${1}
        TITLE="${2}"
        ZOOM_RANGE_X_START=${3}
        ZOOM_RANGE_X_END=${4}
        ZOOM_RANGE_Y_START=${5}
        ZOOM_RANGE_Y_END=${6}
        DRAGEN_RTG_ROC=${7}
        BWAMEM_RTG_ROC=${8}
        VG_1000GP_DEFAULT_RTG_ROC=${9}
        VG_1000GP_TRAINED_RTG_ROC=${10}
        VG_PARENTAL_DEFAULT_RTG_ROC=${11}
        VG_PARENTAL_TRAINED_RTG_ROC=${12}
        docker run \
        -e OUT_NAME=${OUT_NAME} \
        -e TITLE="${TITLE}" \
        -e ZOOM_RANGE_X_START=${ZOOM_RANGE_X_START} \
        -e ZOOM_RANGE_X_END=${ZOOM_RANGE_X_END} \
        -e ZOOM_RANGE_Y_START=${ZOOM_RANGE_Y_START} \
        -e ZOOM_RANGE_Y_END=${ZOOM_RANGE_Y_END} \
        -e DRAGEN_RTG_ROC=${DRAGEN_RTG_ROC} \
        -e BWAMEM_RTG_ROC=${BWAMEM_RTG_ROC} \
        -e VG_1000GP_DEFAULT_RTG_ROC=${VG_1000GP_DEFAULT_RTG_ROC} \
        -e VG_1000GP_TRAINED_RTG_ROC=${VG_1000GP_TRAINED_RTG_ROC} \
        -e VG_PARENTAL_DEFAULT_RTG_ROC=${VG_PARENTAL_DEFAULT_RTG_ROC} \
        -e VG_PARENTAL_TRAINED_RTG_ROC=${VG_PARENTAL_TRAINED_RTG_ROC} \
        -v ${PWD}:${HOME} -w ${HOME} realtimegenomics/rtg-tools:3.8.4 \
            rtg rocplot \
            --title "${TITLE}" \
            --plain \
            --zoom=${ZOOM_RANGE_X_START},${ZOOM_RANGE_Y_START},${ZOOM_RANGE_X_END},${ZOOM_RANGE_Y_END} \
            --curve=${BWAMEM_RTG_ROC}=bwamem \
            --curve=${DRAGEN_RTG_ROC}=dragen \
            --curve=${VG_1000GP_DEFAULT_RTG_ROC}=vg_1000gp_default \
            --curve=${VG_1000GP_TRAINED_RTG_ROC}=vg_1000gp_trained \
            --curve=${VG_PARENTAL_DEFAULT_RTG_ROC}=vg_parent_default \
            --curve=${VG_PARENTAL_TRAINED_RTG_ROC}=vg_parent_trained \
            --svg=${OUT_NAME}.svg
    fi
}

WORKDIR=${HOME}/run_genotyping


# Where should temp files go?
mkdir -p "${WORKDIR}"
export TMPDIR="${WORKDIR}/tmp"
mkdir -p "${TMPDIR}"


# Generate the ROC curve plots for High confident region and sample-specific high confident regions based on default-model deeptrio calls.
cd $WORKDIR
for CHILD_NAME in "HG002" "HG005"; do
    for SNP1KG_EXCLUDED in "" "_NOSNP1KG"; do
        if [[ ${SNP1KG_EXCLUDED} == *"_NOSNP1KG"* ]]; then
            if [[ ${CHILD_NAME} == *"HG002"* ]]; then
                ZOOM_RANGE_X_START=0
                ZOOM_RANGE_Y_START=329000
                ZOOM_RANGE_X_END=7000
                ZOOM_RANGE_Y_END=336000
            elif [[ ${CHILD_NAME} == *"HG005"* ]]; then
                ZOOM_RANGE_X_START=0
                ZOOM_RANGE_Y_START=286000
                ZOOM_RANGE_X_END=7000
                ZOOM_RANGE_Y_END=293000
            fi
            TITLE="${CHILD_NAME} v4.2.1 WGS High Confidence non-1000GP Regions DeepTrio Calls"
        else
            if [[ ${CHILD_NAME} == *"HG002"* ]]; then
                ZOOM_RANGE_X_START=0
                ZOOM_RANGE_Y_START=3867000
                ZOOM_RANGE_X_END=8000
                ZOOM_RANGE_Y_END=3875000
            elif [[ ${CHILD_NAME} == *"HG005"* ]]; then
                ZOOM_RANGE_X_START=0
                ZOOM_RANGE_Y_START=3667500
                ZOOM_RANGE_X_END=9000
                ZOOM_RANGE_Y_END=3676500
            fi
            TITLE="${CHILD_NAME} v4.2.1 WGS High Confidence Regions DeepTrio Calls"
        fi
        run_default_rtgrocplot \
            ${CHILD_NAME}${SNP1KG_EXCLUDED}_deeptrio_default \
            ${TITLE} \
            ${ZOOM_RANGE_X_START} ${ZOOM_RANGE_X_END} ${ZOOM_RANGE_Y_START} ${ZOOM_RANGE_Y_END} \
            "rtg_vcfeval_output_${CHILD_NAME}_allhighconfregions${SNP1KG_EXCLUDED}_DRAGEN_DEEPTRIO_DEFAULT/weighted_roc.tsv.gz" \
            "rtg_vcfeval_output_${CHILD_NAME}_allhighconfregions${SNP1KG_EXCLUDED}_BWAMEM_DEEPTRIO_DEFAULT/weighted_roc.tsv.gz" \
            "rtg_vcfeval_output_${CHILD_NAME}_allhighconfregions${SNP1KG_EXCLUDED}_VG_1000GP_DEEPTRIO_DEFAULT/weighted_roc.tsv.gz" \
            "rtg_vcfeval_output_${CHILD_NAME}_allhighconfregions${SNP1KG_EXCLUDED}_VG_PARENTAL_DEEPTRIO_DEFAULT/weighted_roc.tsv.gz" \
    done
done

# Generate the ROC curve plots for High confident region and sample-specific high confident regions based on trained-model deeptrio calls and evaluated on chromosome 20.
cd $WORKDIR
for CHILD_NAME in "HG002" "HG005"; do
    for SNP1KG_EXCLUDED in "" "_NOSNP1KG"; do
        if [[ ${SNP1KG_EXCLUDED} == *"_NOSNP1KG"* ]]; then
            if [[ ${CHILD_NAME} == *"HG002"* ]]; then
                ZOOM_RANGE_X_START=0
                ZOOM_RANGE_Y_START=5910
                ZOOM_RANGE_X_END=150
                ZOOM_RANGE_Y_END=6060
            elif [[ ${CHILD_NAME} == *"HG005"* ]]; then
                ZOOM_RANGE_X_START=0
                ZOOM_RANGE_Y_START=4300
                ZOOM_RANGE_X_END=200
                ZOOM_RANGE_Y_END=4500
            fi
            TITLE="${CHILD_NAME} v4.2.1 Chromosome 20 High Confidence non-1000GP Regions Trained DeepTrio Calls"
        else
            if [[ ${CHILD_NAME} == *"HG002"* ]]; then
                ZOOM_RANGE_X_START=0
                ZOOM_RANGE_Y_START=82220
                ZOOM_RANGE_X_END=180
                ZOOM_RANGE_Y_END=82400
            elif [[ ${CHILD_NAME} == *"HG005"* ]]; then
                ZOOM_RANGE_X_START=0
                ZOOM_RANGE_Y_START=76700
                ZOOM_RANGE_X_END=250
                ZOOM_RANGE_Y_END=76950
            fi
            TITLE="${CHILD_NAME} v4.2.1 Chromosome 20 High Confidence Regions Trained DeepTrio Calls"
        fi
        run_trained_chr20_rtgrocplot \
            ${CHILD_NAME}${SNP1KG_EXCLUDED}_deeptrio_trained_chr20 \
            ${TITLE} \
            ${ZOOM_RANGE_X_START} ${ZOOM_RANGE_X_END} ${ZOOM_RANGE_Y_START} ${ZOOM_RANGE_Y_END} \
            "rtg_vcfeval_output_${CHILD_NAME}_allhighconfregions${SNP1KG_EXCLUDED}_DRAGEN_DEEPTRIO_DEFAULT_chr20/weighted_roc.tsv.gz" \
            "rtg_vcfeval_output_${CHILD_NAME}_allhighconfregions${SNP1KG_EXCLUDED}_BWAMEM_DEEPTRIO_DEFAULT_chr20/weighted_roc.tsv.gz" \
            "rtg_vcfeval_output_${CHILD_NAME}_allhighconfregions${SNP1KG_EXCLUDED}_VG_1000GP_DEEPTRIO_DEFAULT_chr20/weighted_roc.tsv.gz" \
            "rtg_vcfeval_output_${CHILD_NAME}_allhighconfregions${SNP1KG_EXCLUDED}_VG_1000GP_DEEPTRIO_TRAINED_chr20/weighted_roc.tsv.gz" \
            "rtg_vcfeval_output_${CHILD_NAME}_allhighconfregions${SNP1KG_EXCLUDED}_VG_PARENTAL_DEEPTRIO_DEFAULT_chr20/weighted_roc.tsv.gz" \
            "rtg_vcfeval_output_${CHILD_NAME}_allhighconfregions${SNP1KG_EXCLUDED}_VG_PARENTAL_DEEPTRIO_TRAINED_chr20/weighted_roc.tsv.gz" \
    done
done
