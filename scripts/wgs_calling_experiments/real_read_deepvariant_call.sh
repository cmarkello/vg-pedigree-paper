#!/usr/bin/env bash
# real_read_deeptrio_call.sh: indel realign input BAMs and run the DeepTrio Genotyper

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

INDEL_REALIGN_BAMS=true
WORKDIR=${HOME}/run_genotyping
OUTNAME="${1}"
CHILD_NAME="${2}"
BAM_FILE_CHILD="${3}"
INDEL_REALIGN_BAMS=${4}
INPUT_CHILD_BAM="${CHILD_NAME}.bam"
REF_FASTA="GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.gz"
SEQ_DICT="GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.dict"


# Where should temp files go?
mkdir -p "${WORKDIR}"
export TMPDIR="${WORKDIR}/tmp"
mkdir -p "${TMPDIR}"

# Download data input data
cd $WORKDIR
copy ${BAM_FILE_CHILD} "${WORKDIR}/${INPUT_CHILD_BAM}"
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.gz "${WORKDIR}/${REF_FASTA}"
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.dict -O "${WORKDIR}/${SEQ_DICT}"

# Sort and reorder and indel-realign input trio BAMS
if [[ $INDEL_REALIGN_BAMS != true ]]; then 
    mv "${WORKDIR}/${INPUT_CHILD_BAM}" "${WORKDIR}/indel_realigned.${CHILD_NAME}.bam"
else
    cd $WORKDIR
    SAMPLE_NAME="${CHILD_NAME}"
    INPUT_BAM=${SAMPLE_NAME}.bam
    docker run \
    -e INPUT_BAM=${INPUT_BAM} \
    -v ${PWD}:${HOME} -w ${HOME} quay.io/ucsc_cgl/samtools:latest \
        samtools sort -@ 32 ${INPUT_BAM} -O BAM > positionsorted.${INPUT_BAM} && rm ${INPUT_BAM}

    docker run \
    -e INPUT_BAM=${INPUT_BAM} \
    -e SEQ_DICT=${SEQ_DICT} \
    -v ${PWD}:${HOME} -w ${HOME} broadinstitute/picard:2.21.9 \
      java -Xmx20g -XX:ParallelGCThreads=16 -jar /usr/picard/picard.jar \
      ReorderSam \
      VALIDATION_STRINGENCY=SILENT \
      INPUT=positionsorted.${INPUT_BAM} \
      OUTPUT=reordered.positionsorted.${INPUT_BAM} \
      SEQUENCE_DICTIONARY=${SEQ_DICT} && rm positionsorted.${INPUT_BAM}

    # Indel realign input BAM
    cd $WORKDIR
    docker run \
    -e INPUT_BAM=${INPUT_BAM} \
    -e SAMPLE_NAME=${SAMPLE_NAME} \
    -v ${PWD}:${HOME} -w ${HOME} quay.io/ucsc_cgl/samtools:latest \
    samtools addreplacerg \
    -@ 32 -O BAM -r ID:1 -r LB:lib1 -r SM:${SAMPLE_NAME} -r PL:illumina -r PU:unit1 \
    reordered.positionsorted.${INPUT_BAM} > gatk_ready.reordered.positionsorted.${INPUT_BAM}

    docker run \
    -e INPUT_BAM=${INPUT_BAM} \
    -v ${PWD}:${HOME} -w ${HOME} quay.io/ucsc_cgl/samtools:latest \
    samtools index -@ 32 gatk_ready.reordered.positionsorted.${INPUT_BAM}

    docker run \
    -e SAMPLE_NAME=${SAMPLE_NAME} \
    -e REF_FASTA=${REF_FASTA} \
    -e INPUT_BAM=${INPUT_BAM} \
    -v ${PWD}:${HOME} -w ${HOME} broadinstitute/gatk3:3.8-1 \
      java -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator \
      -R ${REF_FASTA} \
      -I gatk_ready.reordered.positionsorted.${INPUT_BAM} -o ${SAMPLE_NAME}.intervals

    awk -F '[:-]' 'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}' ${SAMPLE_NAME}.intervals > ${SAMPLE_NAME}.intervals.bed
    docker run \
    -e SAMPLE_NAME=${SAMPLE_NAME} \
    -e INPUT_BAM=${INPUT_BAM} \
    -e REF_FASTA=${REF_FASTA} \
    -v ${PWD}:${HOME} -w ${HOME} quay.io/biocontainers/abra2:2.24--h7d875b9_0 \
      --targets ${SAMPLE_NAME}.intervals.bed \
      --in gatk_ready.reordered.positionsorted.${INPUT_BAM} \
      --out indel_realigned.${INPUT_BAM} \
      --ref ${REF_FASTA} \
      --threads 32
    
    docker run \
    -e SAMPLE_NAME=${SAMPLE_NAME} \
    -v ${PWD}:${HOME} -w ${HOME} quay.io/ucsc_cgl/samtools:latest \
        samtools index -@ 32 indel_realigned.${SAMPLE_NAME}.bam
fi


# Run DeepTrio genotyper on input trio BAMs
cd ${WORKDIR}
mkdir -p ${WORKDIR}/tmp_deepvariant
docker run \
-e REF_FASTA=${REF_FASTA} \
-e CHILD_NAME=${CHILD_NAME} \
-e OUTNAME=${OUTNAME} \
-v ${PWD}:${HOME} -w ${HOME} docker://google/deepvariant:1.1.0 \
/bin/bash -c '/opt/deepvariant/bin/run_deepvariant \
--make_examples_extra_args "min_mapping_quality=1" \
--model_type=WGS \
--ref=${REF_FASTA} \
--reads=indel_realigned.${CHILD_NAME}.bam \
--output_vcf=${CHILD_NAME}_DEEPVARIANT.abra_gatk_targets.${OUTNAME}.vcf.gz \
--output_gvcf=${CHILD_NAME}_DEEPVARIANT.abra_gatk_targets.${OUTNAME}.g.vcf.gz \
--intermediate_results_dir=tmp_deepvariant.indel_realigned \
--num_shards=32'

