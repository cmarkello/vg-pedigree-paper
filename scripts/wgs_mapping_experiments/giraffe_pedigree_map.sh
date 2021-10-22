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
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/liftover_snp1kg_grch38_nosegdup.xg -O ${WORKFLOW_INPUT_DIR}/snp1kg_decoys.xg
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/liftover_snp1kg_grch38_nosegdup.gbwt -O ${WORKFLOW_INPUT_DIR}/snp1kg_decoys.gbwt
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/liftover_snp1kg_grch38_nosegdup.gg -O ${WORKFLOW_INPUT_DIR}/snp1kg_decoys.gg
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/liftover_snp1kg_grch38_nosegdup.min -O ${WORKFLOW_INPUT_DIR}/snp1kg_decoys.min
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/liftover_snp1kg_grch38_nosegdup.dist -O ${WORKFLOW_INPUT_DIR}/snp1kg_decoys.dist
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.gz -O ${WORKFLOW_INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.gz
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.fai -O ${WORKFLOW_INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.fai
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.dict -O ${WORKFLOW_INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.dict
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/snpEff_v5_0_GRCh38.99.zip -O ${WORKFLOW_INPUT_DIR}/snpEff_v5_0_GRCh38.99.zip
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/eagle_data_grch38.tar.gz -O ${WORKFLOW_INPUT_DIR}/eagle_data_grch38.tar.gz
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/dt-giraffe-child-0806.tar.gz -O ${WORKFLOW_INPUT_DIR}/dt-giraffe-child-0806.tar.gz
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/dt-giraffe-parent-0806.tar.gz -O ${WORKFLOW_INPUT_DIR}/dt-giraffe-parent-0806.tar.gz
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/dv-giraffe-0507.tar.gz -O ${WORKFLOW_INPUT_DIR}/dv-giraffe-0507.tar.gz
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/HG002.ped -O ${WORKFLOW_INPUT_DIR}/HG002.ped
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/HG005.ped -O ${WORKFLOW_INPUT_DIR}/HG005.ped

INPUT_CHILD_MODEL=""
INPUT_PARENT_MODEL=""
INPUT_DEEPVARIANT_MODEL=""
if [[ $MODEL_USED == *"trained"* ]]; then
    wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/dt-giraffe-child-0806.tar.gz -O ${WORKFLOW_INPUT_DIR}/dt-giraffe-child-0806.tar.gz
    wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/dt-giraffe-parent-0806.tar.gz -O ${WORKFLOW_INPUT_DIR}/dt-giraffe-parent-0806.tar.gz
    wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/dv-giraffe-0507.tar.gz -O ${WORKFLOW_INPUT_DIR}/dv-giraffe-0507.tar.gz
    INPUT_CHILD_MODEL="--deeptrio_child_model ${WORKFLOW_INPUT_DIR}/dt-giraffe-child-0806.tar.gz"
    INPUT_PARENT_MODEL="--deeptrio_parent_model ${WORKFLOW_INPUT_DIR}/dt-giraffe-parent-0806.tar.gz"
    INPUT_DEEPVARIANT_MODEL="--deepvariant_model ${WORKFLOW_INPUT_DIR}/dv-giraffe-0507.tar.gz"
fi

# Download input reads
CHILD_ID=""
MATERNAL_ID=""
PATERNAL_ID=""
if [[ ${CHILD_NAME} == *"HG002"* ]]; then
    CHILD_ID="HG002"
    MATERNAL_ID="HG004"
    PATERNAL_ID="HG003"
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
    wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG005.novaseq.pcr-free.30x.R1.fastq.gz ${WORKFLOW_INPUT_DIR}/${CHILD_ID}.R1.fastq.gz
    wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG005.novaseq.pcr-free.30x.R2.fastq.gz ${WORKFLOW_INPUT_DIR}/${CHILD_ID}.R2.fastq.gz
    wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG006.novaseq.pcr-free.30x.R1.fastq.gz ${WORKFLOW_INPUT_DIR}/${PATERNAL_ID}.R1.fastq.gz
    wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG006.novaseq.pcr-free.30x.R2.fastq.gz ${WORKFLOW_INPUT_DIR}/${PATERNAL_ID}.R2.fastq.gz
    wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG007.novaseq.pcr-free.30x.R1.fastq.gz ${WORKFLOW_INPUT_DIR}/${MATERNAL_ID}.R1.fastq.gz
    wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG007.novaseq.pcr-free.30x.R2.fastq.gz ${WORKFLOW_INPUT_DIR}/${MATERNAL_ID}.R2.fastq.gz
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
# Run the VG Pedigree workflow for the HG002 cohort
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
--xg_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.xg \
--gbwt_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.gbwt \
--graph_gbwt_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.gg \
--minimizer_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.min \
--distance_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.dist \
--id_ranges ${WORKFLOW_INPUT_DIR}/path_list_whole_genome.txt \
--path_list ${WORKFLOW_INPUT_DIR}/path_list_whole_genome.txt \
--ped_file ${WORKFLOW_INPUT_DIR}/HG002.ped \
--eagle_data ${WORKFLOW_INPUT_DIR}/eagle_data_grch38.tar.gz \
--snpeff_database ${WORKFLOW_INPUT_DIR}/snpEff_v5_0_GRCh38.99.zip \
--bam_output \
--use_decoys \
--indel_realign_bams \
--snpeff_annotation

