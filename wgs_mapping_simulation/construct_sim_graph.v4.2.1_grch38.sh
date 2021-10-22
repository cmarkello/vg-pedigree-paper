#!/usr/bin/env bash
# construct_sim_graph.v4.2.1_grch38.sh: HG002 sample graph construction workflow for read simulation experiments

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

WORKDIR=${HOME}/run_sim_reads
WORKFLOW_INPUT_DIR=${WORKDIR}/inputs

# Where should temp files go?
mkdir -p "${WORKDIR}" "${WORKFLOW_INPUT_DIR}"
export TMPDIR="${WORKDIR}/tmp"
mkdir -p "${TMPDIR}"

# Download data input data
cd $WORKDIR
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.gz -O ${WORKFLOW_INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.gz
wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/SupplementaryFiles/HG002_GRCh38_1_22_v4.2.1_benchmark_phased_MHCassembly_StrandSeqANDTrio.vcf.gz -O ${WORKFLOW_INPUT_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark_phased_MHCassembly_StrandSeqANDTrio.vcf.gz

cd ${WORKFLOW_INPUT_DIR}
docker run \
-v ${PWD}:${HOME} -w ${HOME} quay.io/biocontainers/vt:0.57721--heae7c10_3 \
vt view \
-f "ALT~~'\*'" \
HG002_GRCh38_1_22_v4.2.1_benchmark_phased_MHCassembly_StrandSeqANDTrio.vcf.gz \
-o HG002_GRCh38_1_22_v4.2.1_benchmark_phased_MHCassembly_StrandSeqANDTrio.no_asterisk.vcf.gz

# Setup toil-vg virtual environment
cd $WORKDIR
git clone --single-branch --branch vg_pedigree_workflow_deepvariant_grch38 https://github.com/vgteam/toil-vg.git
git clone https://github.com/cmarkello/toil.git
python3 -m venv toilvg_venv
virtualenv --system-site-packages --python python3 venv
source venv/bin/activate
pip3 install ./toil
pip3 install ./toil-vg
toil-vg generate-config --whole_genome >config.cfg
sed -i'' config.cfg -e "s/gcsa-index-mem: '110G'/gcsa-index-mem: '700G'/g" -e "s/gcsa-index-disk: '2200G'/gcsa-index-disk: '4096G'/g" -e "s/construct-mem: '64G'/construct-mem: '128G'/g" -e "s/construct-disk: '64G'/construct-disk: '128G'/g" -e "s/gbwt-index-mem: '35G'/gbwt-index-mem: '70G'/g" -e "s/gbwt-index-disk: '100G'/gbwt-index-disk: '200G'/g"

# Setup construction workflow
export XDG_RUNTIME_DIR=""
TOILVG_OUTSTORE="${WORKDIR}/vg-construct-sim-v4.2.1_grch38-outstore"
TOILVG_JOBSTORE="${WORKDIR}/vg-construct-sim-v4.2.1_grch38-jobstore"
LOGFILE="${WORKDIR}/construct_sim_baseline.v4.2.1_grch38.log"
mkdir -p ${TOILVG_OUTSTORE}
rm -f ${LOGFILE}
export TOIL_SLURM_ARGS="-t 40:00:00"
export SINGULARITY_CACHEDIR=${WORKDIR}/singularity-cache

# Run the construction workflow
toil clean ${TOILVG_JOBSTORE}
toil-vg construct \
${TOILVG_JOBSTORE} \
${TOILVG_OUTSTORE} \
--batchSystem singleMachine \
--container Docker \
--vg_docker quay.io/vgteam/vg:ci-2890-655a9622c3d60e87f14b88d943fbd8554214a975 \
--realTimeLogging \
--logInfo \
--workDir ${WORK_DIR}/tmp \
--cleanWorkDir onSuccess \
--fasta ${WORKFLOW_INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.gz \
--out_name hg002_sample_grch38 \
--logFile ${LOGFILE} \
--xg_index --gbwt_index \
--force_phasing True \
--pangenome \
--sample_graph HG002 \
--merge_graphs \
--keep_vcfs \
--fasta_regions \
--vcf ${WORKFLOW_INPUT_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark_phased_MHCassembly_StrandSeqANDTrio.no_asterisk.vcf.gz \
--vcf_phasing ${WORKFLOW_INPUT_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark_phased_MHCassembly_StrandSeqANDTrio.no_asterisk.vcf.gz \
--statePollingWait 120 \
--rescueJobsFrequency 120

