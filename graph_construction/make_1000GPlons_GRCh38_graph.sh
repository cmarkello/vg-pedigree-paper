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
# VCF order must be 1-22, X, Y to match FASTA
CONST_VCFS=()
CONST_VCFS=($((seq 1 22) | xargs -I {} echo ${WORKFLOW_INPUT_DIR}/ALL.chr{}_GRCh38.genotypes.20170504.no_segdups_gt10kb.vcf.gz))
CONST_VCFS+=(${WORKFLOW_INPUT_DIR}/ALL.chrX_GRCh38.genotypes.20170504.no_segdups_gt10kb.vcf.gz)
CONST_VCFS+=(${WORKFLOW_INPUT_DIR}/ALL.chrY_GRCh38.genotypes.20170504.no_segdups_gt10kb.vcf.gz)
CONST_VCFS_PHASED=()
CONST_VCFS_PHASED=($((seq 1 22) | xargs -I {} echo ${WORKFLOW_INPUT_DIR}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{}.filtered.shapeit2-duohmm-phased.vcf.gz))
CONST_VCFS_PHASED+=(${WORKFLOW_INPUT_DIR}/ALL.chrX_GRCh38.genotypes.20170504.no_segdups_gt10kb.vcf.gz)
TOILVG_OUTSTORE="${WORKDIR}/construct-grch38-liftover-nosegdup-outstore"
TOILVG_JOBSTORE="${WORKDIR}/construct-grch38-liftover-nosegdup-jobstore"
LOGFILE="${WORKDIR}/construct_grch38_liftover_nosegdup.log"
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
--workDir ${WORKDIR}/tmp \
--cleanWorkDir onSuccess \
--config config.cfg \
--fasta ${WORKFLOW_INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.gz \
--out_name liftover_snp1kg_grch38_nosegdup \
--logFile ${LOGFILE} \
--xg_index --gbwt_index --trivial_snarls_index \
--force_phasing True \
--merge_graphs \
--fasta_regions \
--pangenome \
--vcf ${CONST_VCFS[@]} \
--vcf_phasing ${CONST_VCFS_PHASED[@]} \
--statePollingWait 120 \
--rescueJobsFrequency 120 \
--defaultMemory 50Gi --defaultDisk 400G

# Build sampled GBWT-based indexes
cd $WORKDIR/construct-grch38-liftover-nosegdup-outstore
docker run \
-v ${PWD}:${HOME} -w ${HOME} quay.io/vgteam/vg:ci-2351-a64b70c1f9345f0821e3f3a600eb8bbf4fe44bf2 \
    vg index \
    -s liftover_snp1kg_grch38_nosegdup.trivial.snarls \
    -j liftover_snp1kg_grch38_nosegdup.trivial.snarls.dist \
    -x liftover_snp1kg_grch38_nosegdup.xg \
    -t 32 -p

# Build sampled gbwt and gbwt graph
docker run \
-v ${PWD}:${HOME} -w ${HOME} quay.io/vgteam/vg:ci-2351-a64b70c1f9345f0821e3f3a600eb8bbf4fe44bf2 \
    vg gbwt \
    -x liftover_snp1kg_grch38_nosegdup.xg \
    -o liftover_snp1kg_grch38_nosegdup.sampled.gbwt \ 
    -g liftover_snp1kg_grch38_nosegdup.sampled.gg \
    liftover_snp1kg_grch38_nosegdup.gbwt \
    -l -n 64 \
    -t 32 -p

# Build the minimizer index
docker run \
-v ${PWD}:${HOME} -w ${HOME} quay.io/vgteam/vg:ci-2351-a64b70c1f9345f0821e3f3a600eb8bbf4fe44bf2 \
    vg minimizer \
    -g liftover_snp1kg_grch38_nosegdup.sampled.gbwt \
    -i liftover_snp1kg_grch38_nosegdup.sampled.trivial_snarls_dist.min \
    -d liftover_snp1kg_grch38_nosegdup.trivial.snarls.dist \
    -G liftover_snp1kg_grch38_nosegdup.sampled.gg \
    -t 32 -p

mv liftover_snp1kg_grch38_nosegdup.gbwt liftover_snp1kg_grch38_nosegdup.old.gbwt
mv liftover_snp1kg_grch38_nosegdup.sampled.gbwt liftover_snp1kg_grch38_nosegdup.gbwt
mv liftover_snp1kg_grch38_nosegdup.sampled.trivial_snarls_dist.min liftover_snp1kg_grch38_nosegdup.min
mv liftover_snp1kg_grch38_nosegdup.sampled.gg liftover_snp1kg_grch38_nosegdup.gg
mv liftover_snp1kg_grch38_nosegdup.trivial.snarls.dist liftover_snp1kg_grch38_nosegdup.dist

