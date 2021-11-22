#!/usr/bin/env bash
# giraffe_pedigree_map.sh: Main VG Pedigree mapping workflow.

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

function copy() {
    if [ ! -e "${2}" ] ; then
        cp "${1}" "${2}"
    fi
}

CHILD_NAME="${1}"
MODEL_USED="${2}" # "default" or "trained"
WORKDIR=${HOME}/run_giraffe_pedigree_mapping
WORKFLOW_INPUT_DIR=${WORKDIR}/inputs

# Where should temp files go?
mkdir -p "${WORKDIR}" "${WORKFLOW_INPUT_DIR}"

# Download input data
cd $WORKDIR
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/path_list_whole_genome.txt -O ${WORKFLOW_INPUT_DIR}/path_list_whole_genome.txt
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.gz -O ${WORKFLOW_INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.gz
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.fai -O ${WORKFLOW_INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.fai
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.dict -O ${WORKFLOW_INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.dict
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/snpEff_v5_0_GRCh38.99.zip -O ${WORKFLOW_INPUT_DIR}/snpEff_v5_0_GRCh38.99.zip
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/liftover_snp1kg_grch38_nosegdup_noceph.xg -O ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup_noceph.xg
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/liftover_snp1kg_grch38_nosegdup_noceph.gbwt -O ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup_noceph.gbwt
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/liftover_snp1kg_grch38_nosegdup_noceph.gg -O ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup_noceph.gg
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/liftover_snp1kg_grch38_nosegdup_noceph.min -O ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup_noceph.min
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/liftover_snp1kg_grch38_nosegdup_noceph.dist -O ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup_noceph.dist
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/liftover_snp1kg_grch38_nosegdup.xg -O ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.xg
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/liftover_snp1kg_grch38_nosegdup.gbwt -O ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.gbwt
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/liftover_snp1kg_grch38_nosegdup.gg -O ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.gg
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/liftover_snp1kg_grch38_nosegdup.min -O ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.min
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/liftover_snp1kg_grch38_nosegdup.dist -O ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.dist
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/eagle_data_grch38_noceph.tar.gz -O ${WORKFLOW_INPUT_DIR}/eagle_data_grch38_noceph.tar.gz
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/eagle_data_grch38.tar.gz -O ${WORKFLOW_INPUT_DIR}/eagle_data_grch38.tar.gz

INPUT_CHILD_MODEL=""
INPUT_PARENT_MODEL=""
INPUT_DEEPVARIANT_MODEL=""
if [[ $MODEL_USED == *"trained"* ]]; then
    if [[ ${CHILD_NAME} == *"HG001"* ]]; then
        wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/dt-giraffe-child-0821.tar.gz -O ${WORKFLOW_INPUT_DIR}/dt-giraffe-child-0821.tar.gz
        wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/dt-giraffe-parent-0821.tar.gz -O ${WORKFLOW_INPUT_DIR}/dt-giraffe-parent-0821.tar.gz
        INPUT_CHILD_MODEL="--deeptrio_child_model ${WORKFLOW_INPUT_DIR}/dt-giraffe-child-0821.tar.gz"
        INPUT_PARENT_MODEL="--deeptrio_parent_model ${WORKFLOW_INPUT_DIR}/dt-giraffe-parent-0821.tar.gz"
    else
        wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/dt-giraffe-child-0806.tar.gz -O ${WORKFLOW_INPUT_DIR}/dt-giraffe-child-0806.tar.gz
        wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/dt-giraffe-parent-0806.tar.gz -O ${WORKFLOW_INPUT_DIR}/dt-giraffe-parent-0806.tar.gz
        INPUT_CHILD_MODEL="--deeptrio_child_model ${WORKFLOW_INPUT_DIR}/dt-giraffe-child-0806.tar.gz"
        INPUT_PARENT_MODEL="--deeptrio_parent_model ${WORKFLOW_INPUT_DIR}/dt-giraffe-parent-0806.tar.gz"
    fi
    wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/dv-giraffe-0507.tar.gz -O ${WORKFLOW_INPUT_DIR}/dv-giraffe-0507.tar.gz
    INPUT_DEEPVARIANT_MODEL="--deepvariant_model ${WORKFLOW_INPUT_DIR}/dv-giraffe-0507.tar.gz"
fi

# Download input reads
CHILD_ID=""
MATERNAL_ID=""
PATERNAL_ID=""
INPUT_PED_FILE=""
INPUT_EAGLE_DATA=""
INPUT_XG_INDEX=""
INPUT_GBWT_INDEX=""
INPUT_GG_INDEX=""
INPUT_MIN_INDEX=""
INPUT_DIST_INDEX=""
if [[ ${CHILD_NAME} == *"HG002"* ]]; then
    CHILD_ID="HG002"
    MATERNAL_ID="HG004"
    PATERNAL_ID="HG003"
    INPUT_PED_FILE="--ped_file ${WORKFLOW_INPUT_DIR}/HG002.ped"
    INPUT_EAGLE_DATA="--eagle_data ${WORKFLOW_INPUT_DIR}/eagle_data_grch38.tar.gz"
    INPUT_XG_INDEX="--xg_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.xg"
    INPUT_GBWT_INDEX="--gbwt_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.gbwt"
    INPUT_GG_INDEX="--graph_gbwt_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.gg"
    INPUT_MIN_INDEX="--minimizer_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.min"
    INPUT_DIST_INDEX="--distance_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.dist"
    wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/HG002.ped -O ${WORKFLOW_INPUT_DIR}/HG002.ped
    wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG002.novaseq.pcr-free.35x.R1.fastq.gz ${WORKFLOW_INPUT_DIR}/${CHILD_ID}.R1.fastq.gz
    wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG002.novaseq.pcr-free.35x.R2.fastq.gz ${WORKFLOW_INPUT_DIR}/${CHILD_ID}.R2.fastq.gz
    wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG003.novaseq.pcr-free.35x.R1.fastq.gz ${WORKFLOW_INPUT_DIR}/${PATERNAL_ID}.R1.fastq.gz
    wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG003.novaseq.pcr-free.35x.R2.fastq.gz ${WORKFLOW_INPUT_DIR}/${PATERNAL_ID}.R2.fastq.gz
    wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG004.novaseq.pcr-free.35x.R1.fastq.gz ${WORKFLOW_INPUT_DIR}/${MATERNAL_ID}.R1.fastq.gz
    wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG004.novaseq.pcr-free.35x.R2.fastq.gz ${WORKFLOW_INPUT_DIR}/${MATERNAL_ID}.R2.fastq.gz
elif [[ ${CHILD_NAME} == *"HG005"* ]]; then
    CHILD_ID="HG005"
    MATERNAL_ID="HG007"
    PATERNAL_ID="HG006"
    INPUT_PED_FILE="--ped_file ${WORKFLOW_INPUT_DIR}/HG005.ped"
    INPUT_EAGLE_DATA="--eagle_data ${WORKFLOW_INPUT_DIR}/eagle_data_grch38.tar.gz"
    INPUT_XG_INDEX="--xg_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.xg"
    INPUT_GBWT_INDEX="--gbwt_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.gbwt"
    INPUT_GG_INDEX="--graph_gbwt_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.gg"
    INPUT_MIN_INDEX="--minimizer_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.min"
    INPUT_DIST_INDEX="--distance_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.dist"
    wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/HG005.ped -O ${WORKFLOW_INPUT_DIR}/HG005.ped
    wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG005.novaseq.pcr-free.30x.R1.fastq.gz ${WORKFLOW_INPUT_DIR}/${CHILD_ID}.R1.fastq.gz
    wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG005.novaseq.pcr-free.30x.R2.fastq.gz ${WORKFLOW_INPUT_DIR}/${CHILD_ID}.R2.fastq.gz
    wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG006.novaseq.pcr-free.30x.R1.fastq.gz ${WORKFLOW_INPUT_DIR}/${PATERNAL_ID}.R1.fastq.gz
    wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG006.novaseq.pcr-free.30x.R2.fastq.gz ${WORKFLOW_INPUT_DIR}/${PATERNAL_ID}.R2.fastq.gz
    wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG007.novaseq.pcr-free.30x.R1.fastq.gz ${WORKFLOW_INPUT_DIR}/${MATERNAL_ID}.R1.fastq.gz
    wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG007.novaseq.pcr-free.30x.R2.fastq.gz ${WORKFLOW_INPUT_DIR}/${MATERNAL_ID}.R2.fastq.gz
elif [[ ${CHILD_NAME} == *"HG001"* ]]; then
    CHILD_ID="HG001"
    MATERNAL_ID="NA12891"
    PATERNAL_ID="NA12892"
    INPUT_PED_FILE="--ped_file ${WORKFLOW_INPUT_DIR}/HG001.ped"
    INPUT_EAGLE_DATA="--eagle_data ${WORKFLOW_INPUT_DIR}/eagle_data_grch38_noceph.tar.gz"
    INPUT_XG_INDEX="--xg_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup_noceph.xg"
    INPUT_GBWT_INDEX="--gbwt_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup_noceph.gbwt"
    INPUT_GG_INDEX="--graph_gbwt_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup_noceph.gg"
    INPUT_MIN_INDEX="--minimizer_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup_noceph.min"
    INPUT_DIST_INDEX="--distance_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup_noceph.dist"
    wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG001.ped -O ${WORKFLOW_INPUT_DIR}/HG001.ped
    wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG001.novaseq.pcr-free.30x.R1.fastq.gz ${WORKFLOW_INPUT_DIR}/${CHILD_ID}.R1.fastq.gz
    wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG001.novaseq.pcr-free.30x.R2.fastq.gz ${WORKFLOW_INPUT_DIR}/${CHILD_ID}.R2.fastq.gz
    wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/NA12891.novaseq.pcr-free.30x.R1.fastq.gz ${WORKFLOW_INPUT_DIR}/${PATERNAL_ID}.R1.fastq.gz
    wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/NA12891.novaseq.pcr-free.30x.R2.fastq.gz ${WORKFLOW_INPUT_DIR}/${PATERNAL_ID}.R2.fastq.gz
    wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/NA12892.novaseq.pcr-free.30x.R1.fastq.gz ${WORKFLOW_INPUT_DIR}/${MATERNAL_ID}.R1.fastq.gz
    wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/NA12892.novaseq.pcr-free.30x.R2.fastq.gz ${WORKFLOW_INPUT_DIR}/${MATERNAL_ID}.R2.fastq.gz
fi

git clone --single-branch --branch vg_pedigree_workflow_deepvariant_grch38 https://github.com/vgteam/toil-vg.git 
git clone https://github.com/cmarkello/toil.git
python3 -m venv toilvg_venv
source toilvg_venv/bin/activate
pip install ./toil
pip install ./toil-vg
deactivate

# Run the VG Pedigree workflow using Default and Trained DeepVariant and DeepTrio models

for MODEL in "default" "trained" ; do 
# Run the VG Pedigree workflow for the cohort
source $WORKDIR/toilvg_venv/bin/activate
OUTSTORE="${WORKDIR}/${CHILD_NAME}_${MODEL_USED}_vg_pedigree_outstore"
JOBSTORE="${WORKDIR}/${CHILD_NAME}_${MODEL_USED}_vg_pedigree_jobstore"
LOGFILE="${WORKDIR}/${CHILD_NAME}_${MODEL_USED}_vg_pedigree.log"
export TMPDIR="${WORKDIR}/tmp_${CHILD_NAME}_${MODEL_USED}_vg_pedigree"
export TOIL_SLURM_ARGS='-t 20:00:00'
export SINGULARITY_CACHEDIR=${WORKDIR}/singularity_cache
rm -fr ${LOGFILE} ${OUTSTORE} ${TMPDIR}
mkdir -p ${OUTSTORE} ${TMPDIR} $SINGULARITY_CACHEDIR
cd ${WORKDIR}

toil clean ${JOBSTORE}

time toil-vg pedigree \
--genome_build "GRCh38" \
--retryCount 0 \
--rotatingLogging \
--setEnv PATH=$PATH \
--disableProgress \
--realTimeLogging \
--batchSystem singleMachine \
--statePollingWait 120 \
--rescueJobsFrequency 120 \
--container Docker \
--logInfo \
--logFile ${LOGFILE} \
--workDir ${TMPDIR} \
--cleanWorkDir onSuccess \
--whole_genome_config \
--vg_docker 'quay.io/vgteam/vg:v1.31.0' \
${JOBSTORE} \
${OUTSTORE} \
${CHILD_ID} \
${MATERNAL_ID} \
${PATERNAL_ID} \
--sibling_genders 0 \
--sibling_affected 0 \
--mapper giraffe \
--caller deepvariant \
${INPUT_CHILD_MODEL} \
${INPUT_PARENT_MODEL} \
${INPUT_DEEPVARIANT_MODEL} \
--fastq_proband ${WORKFLOW_INPUT_DIR}/${CHILD_ID}.R1.fastq.gz ${WORKFLOW_INPUT_DIR}/${CHILD_ID}.R2.fastq.gz \
--fastq_maternal ${WORKFLOW_INPUT_DIR}/${MATERNAL_ID}.R1.fastq.gz ${WORKFLOW_INPUT_DIR}/${MATERNAL_ID}.R2.fastq.gz \
--fastq_paternal ${WORKFLOW_INPUT_DIR}/${PATERNAL_ID}.R1.fastq.gz ${WORKFLOW_INPUT_DIR}/${PATERNAL_ID}.R2.fastq.gz \
--ref_fasta ${WORKFLOW_INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.gz \
--ref_fasta_index ${WORKFLOW_INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.fai \
--ref_fasta_dict ${WORKFLOW_INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.dict \
--use_haplotypes \
${INPUT_XG_INDEX} \
${INPUT_GBWT_INDEX} \
${INPUT_GG_INDEX} \
${INPUT_MIN_INDEX} \
${INPUT_DIST_INDEX} \
--id_ranges ${WORKFLOW_INPUT_DIR}/path_list_whole_genome.txt \
--path_list ${WORKFLOW_INPUT_DIR}/path_list_whole_genome.txt \
${INPUT_PED_FILE} \
${INPUT_EAGLE_DATA} \
--snpeff_database ${WORKFLOW_INPUT_DIR}/snpEff_v5_0_GRCh38.99.zip \
--bam_output \
--use_decoys \
--indel_realign_bams \
--snpeff_annotation

