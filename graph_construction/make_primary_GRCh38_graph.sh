#!/usr/bin/env bash
# make_1000GPlons_GRCh38_graph.sh: Main GRCh38-based graph construction workflow

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

WORKDIR=${HOME}/run_graph_construction
WORKFLOW_INPUT_DIR=${WORKDIR}/inputs

# Where should temp files go?
mkdir -p "${WORKDIR}" "${WORKFLOW_INPUT_DIR}"
export TMPDIR="${WORKDIR}/tmp"
mkdir -p "${TMPDIR}"

# Download data input data
cd $WORKDIR
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.gz -O ${WORKFLOW_INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.gz

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
TOILVG_OUTSTORE="${WORKDIR}/construct-grch38-primary-outstore"
TOILVG_JOBSTORE="${WORKDIR}/construct-grch38-primary-jobstore"
LOGFILE="${WORKDIR}/construct_grch38_primary.log"
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
--logInfo \
--vg_docker 'quay.io/vgteam/vg:ci-2351-a64b70c1f9345f0821e3f3a600eb8bbf4fe44bf2' \
--workDir ${WORK_DIR}/tmp_primary \
--cleanWorkDir onSuccess \
--config config.cfg \
--fasta ${WORKFLOW_INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.gz \
--out_name primary_grch38 \
--logFile construct_grch38_fast.log \
--xg_index --snarls_index --trivial_snarls_index --distance_index --id_ranges_index \
--fasta_regions \
--primary \
--statePollingWait 120 \
--rescueJobsFrequency 120 \
--defaultMemory 50Gi --defaultDisk 400G

# Build path-cover GBWT and Minimizer indexes
cd $WORKDIR/construct-grch38-primary-outstore

# Build path-cover-based gbwt and gbwt graph
docker run \
-v ${PWD}:${HOME} -w ${HOME} quay.io/vgteam/vg:ci-2351-a64b70c1f9345f0821e3f3a600eb8bbf4fe44bf2 \
    vg gbwt \
    -x primary_grch38.xg \
    -o primary_grch38.gbwt \
    -g primary_grch38.gg \
    -P \
    -t 32 -p

# Build the minimizer index
docker run \
-v ${PWD}:${HOME} -w ${HOME} quay.io/vgteam/vg:ci-2351-a64b70c1f9345f0821e3f3a600eb8bbf4fe44bf2 \
    vg minimizer \
    -g primary_grch38.gbwt \
    -i primary_grch38.min \
    -d primary_grch38.dist \
    -G primary_grch38.gg \
    -t 32 -p


