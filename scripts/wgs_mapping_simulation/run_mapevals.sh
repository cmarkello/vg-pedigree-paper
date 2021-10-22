#!/bin/bash
# run_mapevals.sh: Evaluate reads simulated from the HG002 sample dataset and region stratifications.

WORKDIR=${HOME}/run_sim_reads


# Where should temp files go?
mkdir -p "${WORKDIR}"
export TMPDIR="${WORKDIR}/tmp"
mkdir -p "${TMPDIR}"

cd ${WORKDIR}

declare -a SIMULATED_READ_DIR_LIST=( "high_conf_hg002_v4.2.1_regions_1M" "all_difficult_regions_hg002_v4.2.1_regions_1M" "alllowmapandsegdupregions_hg002_v4.2.1_regions_100M" "mhc_hg002_v4.2.1_regions_100M" "cmrg_hg002_v4.2.1_regions_100M" "high_conf_NO1000GP_hg002_v4.2.1_regions_1M" "all_difficult_regions_NO1000GP_hg002_v4.2.1_regions_1M" "alllowmapandsegdupregions_NO1000GP_hg002_v4.2.1_regions_100M" "mhc_hg002_NO1000GP_v4.2.1_regions_100M" "cmrg_hg002_NO1000GP_v4.2.1_regions_100M" )
declare -a SIMULATED_READ_DIR_LIST=( "." "high_conf_hg002_v4.2.1_regions" "all_difficult_regions_hg002_v4.2.1_regions" "alllowmapandsegdupregions_hg002_v4.2.1_regions" "mhc_hg002_v4.2.1_regions" "cmrg_hg002_v4.2.1_regions" "high_conf_NO1000GP_hg002_v4.2.1_regions" )

for SIMULATED_READ_DIR in "${SIMULATED_READ_DIR_LIST[@]}" ; do
    cd ${WORK_DIR}

    # RUN BWAMEM
    cd $WORKDIR
    docker run \
    -v ${WORK_DIR}/${SIMULATED_READ_DIR}:${SIMULATED_READ_DIR} \
    -v ${PWD}:${HOME} -w ${HOME} biocontainers/bwa:v0.7.17_cv1 \
        bwa mem \
        -t 16 \
        GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna \
        ${SIMULATED_READ_DIR}/sim.paired.fq
        > mapped.sam
    
    docker run -v ${PWD}:${HOME} -w ${HOME} quay.io/ucsc_cgl/samtools:latest \
        samtools view -F 2048 -b mapped.sam > mapped.primary.bam
    
    docker run -v ${PWD}:${HOME} -w ${HOME} quay.io/ucsc_cgl/samtools:latest \
        samtools view -f 2048 -b mapped.sam > mapped.secondary.bam

    docker run -e SIMULATED_READ_DIR=${SIMULATED_READ_DIR} -v ${PWD}:${HOME} -w ${HOME} docker://quay.io/vgteam/vg:v1.31.0 \
    /bin/bash -c 'set -eo pipefail ; \
    SIMULATED_READ_DIR=${SIMULATED_READ_DIR} ; \
    time vg inject -x hg002_sample_grch38.xg mapped.primary.bam > mapped.primary.bwa.${SIMULATED_READ_DIR}.gam ; \
    time vg inject -x hg002_sample_grch38.xg mapped.secondary.bam > mapped.secondary.bwa.${SIMULATED_READ_DIR}.gam ; \
    time vg view -aj mapped.primary.bwa.${SIMULATED_READ_DIR}.gam | sed "s/\/1/_1/g" | sed "s/\/2/_2/g" | vg view -aGJ - | vg annotate -m -x hg002_sample_grch38.xg -a - | vg gamcompare -r 100 -s - ${SIMULATED_READ_DIR}/sim.gam 2> count | vg view -aj - > compared.primary.json ; \
    time vg view -aj mapped.secondary.bwa.${SIMULATED_READ_DIR}.gam | sed "s/\/1/_1/g" | sed "s/\/2/_2/g" | vg view -aGJ - | vg annotate -m -x hg002_sample_grch38.xg -a - | vg gamcompare -r 100 - ${SIMULATED_READ_DIR}/sim.gam | vg view -aj - > compared.secondary.json'

    python ${HOME}/vg_pedigree_manuscript_methods/wgs_mapping_simulation/combine_reads.py compared.primary.json compared.secondary.json bwamem.compared.json

    docker run -e SIMULATED_READ_DIR=${SIMULATED_READ_DIR} -v ${PWD}:${HOME} -w ${HOME} docker://quay.io/vgteam/vg:v1.31.0 \
    /bin/bash -c 'set -o pipefail ; \
    sed -i "/^$/d" bwamem.compared.json ; \
    printf "correct\tmq\tscore\taligner\n" > bwamem_roc_stats.${SIMULATED_READ_DIR}.tsv ; \
    READS="HG002" ; \
    PAIRING="paired" ; \
    GRAPH="hg002_sample_grch38" ; \
    ALGORITHM="bwamem" ; \
    CORRECT_COUNT="$(grep correctly_mapped bwamem.compared.json | wc -l)"; \
    SCORE="$(sed -n "2p" count | sed "s/[^0-9\.]//g")"; \
    MAPQ="$(grep mapping_quality\":\ 60 bwamem.compared.json | wc -l)"; \
    echo test ; \
    MAPQ60="$(grep -v correctly_mapped bwamem.compared.json | grep mapping_quality\":\ 60 | wc -l)"; \
    echo test2 ; \
    IDENTITY=$(jq ".identity" bwamem.compared.json | awk "{sum+=\$1} END {print sum/NR}"); \
    echo ${GRAPH} ${READS} ${PAIRING} ${CORRECT_COUNT} ${MAPQ} ${MAPQ60} ${SCORE}; \
    printf "graph\talgorithm\treads\tspeed\tpairing\tcorrect\tmapq60\twrong_mapq60\tidentity\tscore\n" > report_bwamem.${SIMULATED_READ_DIR}.tsv; \
    printf "${GRAPH}\t${ALGORITHM}\t${READS}\t${PAIRING}\t-\t${CORRECT_COUNT}\t${MAPQ}\t${MAPQ60}\t${IDENTITY}\t${SCORE}\n" >> report_bwamem.${SIMULATED_READ_DIR}.tsv; \
    jq -r "(if .correctly_mapped then 1 else 0 end|tostring) + \",\" + (.mapping_quality|tostring) + \",\" + (.score|tostring)" bwamem.compared.json | sed "s/\,/\t/g" | sed "s/$/\tbwamem_${GRAPH}${READS}${PAIRING}/" | sed "s/single//g ; s/paired/-pe/g ; s/null/0/g" >> bwamem_roc_stats.${SIMULATED_READ_DIR}.tsv'

    # RUN GIRAFFE PRIMARY
    docker run -e SIMULATED_READ_DIR=${SIMULATED_READ_DIR} -v ${HOME}/run_graph_construction/construct-grch38-primary-outstore:${HOME}/construct-grch38-primary-outstore -v ${PWD}:${HOME} -w ${HOME} docker://quay.io/vgteam/vg:v1.31.0 \
    /bin/bash -c 'set -o pipefail ; \
    PARAM_PRESET="default" ; \
    PAIRED="-i" ; \
    time vg giraffe -x construct-grch38-primary-outstore/primary_grch38.xg -H construct-grch38-primary-outstore/primary_grch38.gbwt -m construct-grch38-primary-outstore/primary_grch38.min -d construct-grch38-primary-outstore/primary_grch38.dist -G ${SIMULATED_READ_DIR}/sim.gam -b ${PARAM_PRESET} ${PAIRED} -t 22 > mapped.giraffe_primary.${SIMULATED_READ_DIR}.gam ; \
    time vg gamcompare -r 100 -s <(vg annotate -m -x construct-grch38-primary-outstore/primary_grch38.xg -a mapped.giraffe_primary.${SIMULATED_READ_DIR}.gam) ${SIMULATED_READ_DIR}/sim.gam 2>count | vg view -aj - > giraffe_primary.compared.json ; \
    GRAPH="hg002_sample_grch38" ; \
    GBWT="sampled.64" ; \
    READS="HG002" ; \
    PARAM_PRESET="default" ; \
    PAIRING="paired" ; \
    CORRECT_COUNT="$(sed -n "1p" count | sed "s/[^0-9]//g")" ; \
    SCORE="$(sed -n "2p" count | sed "s/[^0-9\.]//g")" ; \
    MAPQ="$(grep mapping_quality\":\ 60 giraffe_primary.compared.json | wc -l)" ; \
    MAPQ60="$(grep -v correctly_mapped giraffe_primary.compared.json | grep mapping_quality\":\ 60 | wc -l)" ; \
    IDENTITY=$(jq ".identity" giraffe_primary.compared.json | awk "{sum+=\$1} END {print sum/NR}") ; \
    echo ${GRAPH} ${GBWT} ${READS} ${PARAM_PRESET}${PAIRING} ${CORRECT_COUNT} ${MAPQ} ${MAPQ60} ${IDENTITY} ${SCORE} ; \
    printf "graph\tgbwt\treads\tpairing\tspeed\tcorrect\tmapq60\twrong_mapq60\tidentity\tscore\n" > giraffe_primary_report.${SIMULATED_READ_DIR}.tsv ; \
    printf "${GRAPH}\t${GBWT}\t${READS}\t${PARAM_PRESET}\t${PAIRING}\t${CORRECT_COUNT}\t${MAPQ}\t${MAPQ60}\t${IDENTITY}\t${SCORE}\n" >> giraffe_primary_report.${SIMULATED_READ_DIR}.tsv ; \
    printf "correct\tmq\tscore\taligner\n" > giraffe_primary_roc_stats.${SIMULATED_READ_DIR}.tsv ; \
    jq -r "(if .correctly_mapped then 1 else 0 end|tostring) + \",\" + (.mapping_quality|tostring) + \",\" + (.score|tostring)" giraffe_primary.compared.json | sed "s/\,/\t/g" | sed "s/$/\tgiraffe_primary_${PARAM_PRESET}_${GRAPH}${GBWT}${READS}${PAIRING}/" >> giraffe_primary_roc_stats.${SIMULATED_READ_DIR}.tsv ; \
    sed -i "s/single//g ; s/paired/-pe/g ; s/null/0/g" giraffe_primary_roc_stats.${SIMULATED_READ_DIR}.tsv'

    # RUN GIRAFFE 1000GP
    docker run -e SIMULATED_READ_DIR=${SIMULATED_READ_DIR} -v ${HOME}/run_graph_construction/construct-grch38-liftover-nosegdup-outstore:${HOME}/construct-grch38-liftover-nosegdup-outstore -v ${PWD}:${HOME} -w ${HOME} docker://quay.io/vgteam/vg:v1.31.0 \
    /bin/bash -c 'set -o pipefail ; \
    PARAM_PRESET="default" ; \
    PAIRED="-i" ; \
    time vg giraffe -x construct-grch38-liftover-nosegdup-outstore/liftover_snp1kg_grch38_nosegdup.xg -H construct-grch38-liftover-nosegdup-outstore/liftover_snp1kg_grch38_nosegdup.gbwt -m construct-grch38-liftover-nosegdup-outstore/liftover_snp1kg_grch38_nosegdup.min -d construct-grch38-liftover-nosegdup-outstore/liftover_snp1kg_grch38_nosegdup.dist -G ${SIMULATED_READ_DIR}/sim.gam -b ${PARAM_PRESET} ${PAIRED} -t 22 > mapped.giraffe_snp1kg.${SIMULATED_READ_DIR}.gam ; \
    time vg gamcompare -r 100 -s <(vg annotate -m -x construct-grch38-liftover-nosegdup-outstore/liftover_snp1kg_grch38_nosegdup.xg -a mapped.giraffe_snp1kg.${SIMULATED_READ_DIR}.gam) ${SIMULATED_READ_DIR}/sim.gam 2>count | vg view -aj - > giraffe_snp1kg.compared.json ; \
    GRAPH="hg002_sample_grch38" ; \
    GBWT="sampled.64" ; \
    READS="HG002" ; \
    PARAM_PRESET="default" ; \
    PAIRING="paired" ; \
    CORRECT_COUNT="$(sed -n "1p" count | sed "s/[^0-9]//g")" ; \
    SCORE="$(sed -n "2p" count | sed "s/[^0-9\.]//g")" ; \
    MAPQ="$(grep mapping_quality\":\ 60 giraffe_snp1kg.compared.json | wc -l)" ; \
    MAPQ60="$(grep -v correctly_mapped giraffe_snp1kg.compared.json | grep mapping_quality\":\ 60 | wc -l)" ; \
    IDENTITY=$(jq ".identity" giraffe_snp1kg.compared.json | awk "{sum+=\$1} END {print sum/NR}") ; \
    echo ${GRAPH} ${GBWT} ${READS} ${PARAM_PRESET}${PAIRING} ${CORRECT_COUNT} ${MAPQ} ${MAPQ60} ${IDENTITY} ${SCORE} ; \
    printf "graph\tgbwt\treads\tpairing\tspeed\tcorrect\tmapq60\twrong_mapq60\tidentity\tscore\n" > giraffe_snp1kg_report.${SIMULATED_READ_DIR}.tsv ; \
    printf "${GRAPH}\t${GBWT}\t${READS}\t${PARAM_PRESET}\t${PAIRING}\t${CORRECT_COUNT}\t${MAPQ}\t${MAPQ60}\t${IDENTITY}\t${SCORE}\n" >> giraffe_snp1kg_report.${SIMULATED_READ_DIR}.tsv ; \
    printf "correct\tmq\tscore\taligner\n" > giraffe_snp1kg_roc_stats.${SIMULATED_READ_DIR}.tsv ; \
    jq -r "(if .correctly_mapped then 1 else 0 end|tostring) + \",\" + (.mapping_quality|tostring) + \",\" + (.score|tostring)" giraffe_snp1kg.compared.json | sed "s/\,/\t/g" | sed "s/$/\tgiraffe_snp1kg_${PARAM_PRESET}_${GRAPH}${GBWT}${READS}${PAIRING}/" >> giraffe_snp1kg_roc_stats.${SIMULATED_READ_DIR}.tsv ; \
    sed -i "s/single//g ; s/paired/-pe/g ; s/null/0/g" giraffe_snp1kg_roc_stats.${SIMULATED_READ_DIR}.tsv'

    # RUN GIRAFFE PARENTAL
    docker run -e SIMULATED_READ_DIR=${SIMULATED_READ_DIR} -v ${HOME}/run_giraffe_pedigree_mapping/HG002_default_vg_pedigree_outstore:${HOME}/HG002_vg_pedigree_outstore -v ${PWD}:${HOME} -w ${HOME} docker://quay.io/vgteam/vg:v1.31.0 \
    /bin/bash -c 'set -o pipefail ; \
    PARAM_PRESET="default" ; \
    PAIRED="-i" ; \
    time vg giraffe -x HG002_vg_pedigree_outstore/HG002.parental.graphs.xg -H HG002_vg_pedigree_outstore/HG002.parental.graphs.gbwt -m HG002_vg_pedigree_outstore/HG002.parental.graphs.min -d HG002_vg_pedigree_outstore/HG002.parental.graphs.dist -G ${SIMULATED_READ_DIR}/sim.gam -b ${PARAM_PRESET} ${PAIRED} -t 22 > mapped.giraffe_parental.${SIMULATED_READ_DIR}.gam ; \
    time vg gamcompare -r 100 -s <(vg annotate -m -x HG002_vg_pedigree_outstore/HG002.parental.graphs.xg -a mapped.giraffe_parental.${SIMULATED_READ_DIR}.gam) ${SIMULATED_READ_DIR}/sim.gam 2>count | vg view -aj - > giraffe_parental.compared.json ; \
    GRAPH="hg002_sample_grch38" ; \
    GBWT="sampled.64" ; \
    READS="HG002" ; \
    PARAM_PRESET="default" ; \
    PAIRING="paired" ; \
    CORRECT_COUNT="$(sed -n "1p" count | sed "s/[^0-9]//g")" ; \
    SCORE="$(sed -n "2p" count | sed "s/[^0-9\.]//g")" ; \
    MAPQ="$(grep mapping_quality\":\ 60 giraffe_parental.compared.json | wc -l)" ; \
    MAPQ60="$(grep -v correctly_mapped giraffe_parental.compared.json | grep mapping_quality\":\ 60 | wc -l)" ; \
    IDENTITY=$(jq ".identity" giraffe_parental.compared.json | awk "{sum+=\$1} END {print sum/NR}") ; \
    echo ${GRAPH} ${GBWT} ${READS} ${PARAM_PRESET}${PAIRING} ${CORRECT_COUNT} ${MAPQ} ${MAPQ60} ${IDENTITY} ${SCORE} ; \
    printf "graph\tgbwt\treads\tpairing\tspeed\tcorrect\tmapq60\twrong_mapq60\tidentity\tscore\n" > giraffe_parental_report.${SIMULATED_READ_DIR}.tsv ; \
    printf "${GRAPH}\t${GBWT}\t${READS}\t${PARAM_PRESET}\t${PAIRING}\t${CORRECT_COUNT}\t${MAPQ}\t${MAPQ60}\t${IDENTITY}\t${SCORE}\n" >> giraffe_parental_report.${SIMULATED_READ_DIR}.tsv ; \
    printf "correct\tmq\tscore\taligner\n" > giraffe_parental_roc_stats.${SIMULATED_READ_DIR}.tsv ; \
    jq -r "(if .correctly_mapped then 1 else 0 end|tostring) + \",\" + (.mapping_quality|tostring) + \",\" + (.score|tostring)" giraffe_parental.compared.json | sed "s/\,/\t/g" | sed "s/$/\tgiraffe_parental_${PARAM_PRESET}_${GRAPH}${GBWT}${READS}${PAIRING}/" >> giraffe_parental_roc_stats.${SIMULATED_READ_DIR}.tsv ; \
    sed -i "s/single//g ; s/paired/-pe/g ; s/null/0/g" giraffe_parental_roc_stats.${SIMULATED_READ_DIR}.tsv'
done

