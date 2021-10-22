#!/bin/bash
# sim_reads.HG002.sh: Simulate reads from the HG002 sample graph and stratify by various regions.

set -ex
set -o pipefail

function make_bedfile() {
    if [ ! -e "${3}" ] ; then
        docker run \
        -e HIGH_CONF_BED=${1} \
        -e REGION_BED=${2} \
        -e HIGH_CONF_REGION_BED=${3} \
        -v ${PWD}:${HOME} -w ${HOME} quay.io/biocontainers/bedtools:2.27.0--1
        bedtools intersect \
        -a ${HIGH_CONF_BED} \
        -b ${REGION_BED} \
        > ${HIGH_CONF_REGION_BED}
    fi
}

function make_nosnp1kg_bedfile() {
    if [ ! -e "${3}" ] ; then
        docker run \
        -e HIGH_CONF_BED=${1} \
        -e SNP1KG_VCF_SITES_FILE=${2} \
        -e HIGH_CONF_NOSNP1KG_BED=${3} \
        -v ${PWD}:${HOME} -w ${HOME} quay.io/biocontainers/bedtools:2.27.0--1
        bedtools subtract \
        -a ${HIGH_CONF_BED} \
        -b ${SNP1KG_VCF_SITES_FILE} \
        > ${HIGH_CONF_NOSNP1KG_BED}
    fi
}


WORKDIR=${HOME}/run_sim_reads


# Where should temp files go?
mkdir -p "${WORKDIR}"
export TMPDIR="${WORKDIR}/tmp"
mkdir -p "${TMPDIR}"

# Download data input data
cd $WORKDIR
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG002.novaseq.pcr-free.35x.R1.fastq.gz "${WORKDIR}/${SAMPLE_NAME}.R1.fastq.gz"
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG002.novaseq.pcr-free.35x.R2.fastq.gz "${WORKDIR}/${SAMPLE_NAME}.R2.fastq.gz"
wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed "${WORKDIR}/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/genome-stratifications/v2.0/GRCh38/union/GRCh38_alldifficultregions.bed.gz "${WORKDIR}/GRCh38_alldifficultregions.bed.gz"
wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/genome-stratifications/v2.0/GRCh38/union/GRCh38_alllowmapandsegdupregions.bed.gz "${WORKDIR}/GRCh38_alllowmapandsegdupregions.bed.gz"
wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/genome-stratifications/v2.0/GRCh38/OtherDifficult/GRCh38_MHC.bed.gz "${WORKDIR}/GRCh38_MHC.bed.gz"
wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SmallVariant/HG002_GRCh38_CMRG_smallvar_v1.00.bed "${WORKDIR}/HG002_GRCh38_CMRG_smallvar_v1.00.bed"
wget_download https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz "${WORKDIR}/bedtools-2.30.0.tar.gz" && tar -xzf "${WORKDIR}/bedtools-2.30.0.tar.gz"
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/ALL.GRCh38.genotypes.20170504.rename.sites.sorted.vcf.gz "${WORKDIR}/ALL.GRCh38.genotypes.20170504.rename.sites.sorted.vcf.gz"

# Format the bed files
cd $WORKDIR

make_bedfile HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed ALL.GRCh38.genotypes.20170504.rename.sites.sorted.vcf.gz HG002_snp1kg_liftover_grch38.intersect.bed
docker run -v ${PWD}:${HOME} -w ${HOME} quay.io/biocontainers/bedtools:2.27.0--1 \
/bin/sh -c 'bedtools slop \
-i HG002_snp1kg_liftover_grch38.intersect.bed \
-g bedtools2/genomes/human.hg38.genome \
-b 160 \
> HG002_snp1kg_liftover_grch38.intersect.160bp_slop.bed'

make_bedfile HG002_GRCh38_CMRG_smallvar_v1.00.bed ALL.GRCh38.genotypes.20170504.rename.sites.sorted.vcf.gz HG002_snp1kg_liftover_grch38.CMRG.intersect.bed
docker run -v ${PWD}:${HOME} -w ${HOME} quay.io/biocontainers/bedtools:2.27.0--1 \
/bin/sh -c 'bedtools slop \
-i HG002_snp1kg_liftover_grch38.CMRG.intersect.bed \
-g bedtools2/genomes/human.hg38.genome \
-b 160 \
> HG002_snp1kg_liftover_grch38.CMRG.intersect.160bp_slop.bed'

make_nosnp1kg_bedfile HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed HG002_snp1kg_liftover_grch38.intersect.160bp_slop.bed HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.NO_SNP1KG.bed
make_nosnp1kg_bedfile HG002_GRCh38_CMRG_smallvar_v1.00.bed HG002_snp1kg_liftover_grch38.CMRG.intersect.160bp_slop.bed HG002_GRCh38_CMRG_smallvar_v1.00.NO_SNP1KG.bed

for REGION in "MHC.bed" "alllowmapandsegdupregions.bed" "alldifficultregions.bed" ; do
    make_bedfile HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed GRCh38_${REGION}.gz HG002_GRCh38_v4.2.1.${REGION}
    make_bedfile HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.NO_SNP1KG.bed GRCh38_${REGION}.gz HG002_GRCh38_v4.2.1.NO_SNP1KG.${REGION}
done


# Format the reads
zcat HG002.R1.fastq.gz | seqtk sample -s100 - 1000000 | gzip - > HG002.R1-1m.fq.gz
zcat HG002.R2.fastq.gz | seqtk sample -s100 - 1000000 | gzip - > HG002.R2-1m.fq.gz
paste <(zcat HG002.R1-1m.fq.gz) <(zcat HG002.R2-1m.fq.gz) | paste - - - - | shuf | awk -F'\t' '{OFS="\n"; print $1,$3,$5,$7 > "HG002.R1-shuffled-1m.fq"; print $2,$4,$6,$8 > "HG002.R2-shuffled-1m.fq"}'
gzip HG002.R1-shuffled-1m.fq
gzip HG002.R2-shuffled-1m.fq
docker run -v ${PWD}:${HOME} -w ${HOME} quay.io/biocontainers/bbmap:38.93--he522d1c_0 \
reformat.sh t=2 in=HG002.R1-shuffled-1m.fq.gz in2=HG002.R2-shuffled-1m.fq.gz out=HG002_merged_interleaved-shuffled-1m.fastq.gz

## SIMULATE BASELINE READS

# Simulate 1M reads for all high confident region stratification
docker run \
-e NREADS=1000000 \
-e FASTQ=HG002_merged_interleaved-shuffled-1m.fastq.gz \
-v ${PWD}:${HOME} -w ${HOME} quay.io/vgteam/vg:ci-2890-655a9622c3d60e87f14b88d943fbd8554214a975 \
/bin/bash -c 'vg sim -r -I -n $NREADS -a -s 12345 -p 570 -v 165 -i 0.00029 -x hg002_sample_grch38.xg -g hg002_sample_grch38.gbwt --sample-name HG002 --ploidy-regex "hs38d1:0,chrNC_007605:0,chrX:1,chrY:1,chrY_.*:1,chrEBV:0,.*:2" -F $FASTQ > sim.1m.raw.gam'

# Simulate 100M reads for difficult region stratification
docker run \
-e NREADS=100000000 \
-e FASTQ=HG002_merged_interleaved-shuffled-1m.fastq.gz \
-v ${PWD}:${HOME} -w ${HOME} quay.io/vgteam/vg:ci-2890-655a9622c3d60e87f14b88d943fbd8554214a975 \
/bin/bash -c 'vg sim -r -I -n $NREADS -a -s 12345 -p 570 -v 165 -i 0.00029 -x hg002_sample_grch38.xg -g hg002_sample_grch38.gbwt --sample-name HG002 --ploidy-regex "hs38d1:0,chrNC_007605:0,chrX:1,chrY:1,chrY_.*:1,chrEBV:0,.*:2" -F $FASTQ > sim.100m.raw.gam'

declare -a REGION_LIST=( "high_conf_hg002_v4.2.1_regions_1M" "all_difficult_regions_hg002_v4.2.1_regions_1M" "alllowmapandsegdupregions_hg002_v4.2.1_regions_100M" "mhc_hg002_v4.2.1_regions_100M" "cmrg_hg002_v4.2.1_regions_100M" "high_conf_NO1000GP_hg002_v4.2.1_regions_1M" "all_difficult_regions_NO1000GP_hg002_v4.2.1_regions_1M" "alllowmapandsegdupregions_NO1000GP_hg002_v4.2.1_regions_100M" "mhc_hg002_NO1000GP_v4.2.1_regions_100M" "cmrg_hg002_NO1000GP_v4.2.1_regions_100M" )
declare -a BED_FILE_LIST=( "HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" "HG002_GRCh38_v4.2.1.alldifficultregions.bed" "HG002_GRCh38_v4.2.1.alllowmapandsegdupregions.bed" "HG002_GRCh38_v4.2.1.MHC.bed" "HG002_GRCh38_CMRG_smallvar_v1.00.bed" "HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.NO_SNP1KG.bed" "HG002_GRCh38_v4.2.1.NO_SNP1KG.alldifficultregions.bed" "HG002_GRCh38_v4.2.1.NO_SNP1KG.alllowmapandsegdupregions.bed" "HG002_GRCh38_v4.2.1.NO_SNP1KG.MHC.bed" "HG002_GRCh38_CMRG_smallvar_v1.00.NO_SNP1KG.bed" )

for index in "${!REGION_LIST[@]}"; do
    REGION="${REGION_LIST[index]}"
    BED_FILE="${BED_FILE_LIST[index]}"
    mkdir -p ${WORKDIR}/${REGION}
    cp hg002_sample_grch38.vg ${REGION}/
    if [[ ${REGION} == *"alllowmapandsegdupregions"* || ${REGION} == *"MHC"* || ${REGION} == *"CMRG"* ]]; then
        cp sim.100m.raw.gam ${WORKDIR}/${REGION}/sim.raw.gam
    else
        cp sim.1m.raw.gam ${WORKDIR}/${REGION}/sim.raw.gam
    fi
    cd ${WORKDIR}/${REGION}
    docker run \
    -e REGION=${REGION} \
    -e BED_FILE=${BED_FILE} \
    -v ${WORKDIR}:${HOME}/bed_files \
    -v ${PWD}:${HOME} -w ${HOME} quay.io/vgteam/vg:ci-2890-655a9622c3d60e87f14b88d943fbd8554214a975 \
    /bin/bash -xc 'set -eo pipefail && \
    echo annotate with paths && \
    time vg annotate \
    -p -x hg002_sample_grch38.vg -a sim.raw.gam \
    --bed-name bed_files/${BED_FILE} \
    --threads 32 > sim.${REGION}.gam && \
    time vg filter -i -U -F "" sim.${REGION}.gam > sim.gam && \
    rm -f sim.${REGION}.gam sim.fq.gz sim.fq && \
    echo convert to fastq && \
    time vg view -X -a sim.gam | gzip > sim.fq.gz && \
    gunzip sim.fq.gz && \
    sed "s/_1\$//g" sim.fq | sed "s/_2\$//g" > sim.paired.fq && \
    echo format true position information && \
    time vg view -a sim.gam | jq -c -r "[.name] + if (.annotation.features | length) > 0 then [.annotation.features | join(\",\")] else [\".\"] end + if .refpos != null then [.refpos[] | .name, if .offset != null then .offset else 0 end] else [] end + [.score] + if .mapping_quality == null then [0] else [.mapping_quality] end | @tsv" > true.pos'
done




