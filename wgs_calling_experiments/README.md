# Variant Calling and Evaluation Scripts

These scripts run the DRAGEN and DeepTrio variant callers on the mapped output of BWA-MEM, DRAGEN, VG Giraffe, and VG Pedigree alignments located in `${HOME}/run_bwamem_mapping`, `${HOME}/run_dragen_mapping`, `${HOME}/run_giraffe_mapping`, and `${HOME}/run_giraffe_pedigree_mapping` respectively.

## Running the scripts

First, run the Illumina DRAGEN variant caller on the HG002 and HG005 sample sets of alignments.

```
./real_read_dragen_call.sh HG002 ${HOME}/run_bwamem_mapping/bwamem_HG002.sam false
./real_read_dragen_call.sh HG002 ${HOME}/run_dragen_mapping/dragen_output_HG002.bam false
./real_read_dragen_call.sh HG002 ${HOME}/run_giraffe_mapping/giraffe_HG002.bam true
./real_read_dragen_call.sh HG002 ${HOME}/run_giraffe_pedigree_mapping/HG002_default_vg_pedigree_outstore/HG002_merged.bam true
./real_read_dragen_call.sh HG005 ${HOME}/run_bwamem_mapping/bwamem_HG005.sam false
./real_read_dragen_call.sh HG005 ${HOME}/run_dragen_mapping/dragen_output_HG005.bam false
./real_read_dragen_call.sh HG005 ${HOME}/run_giraffe_mapping/giraffe_HG005.bam true
./real_read_dragen_call.sh HG005 ${HOME}/run_giraffe_pedigree_mapping/HG005_default_vg_pedigree_outstore/HG005_merged.bam true
```

Next, run the DeepTrio variant caller on the HG002 and HG005 trio sample sets of alignments.

```
./real_read_deeptrio_call.sh BWAMEM HG002 HG004 HG003 ${HOME}/run_bwamem_mapping/bwamem_HG002.sam ${HOME}/run_bwamem_mapping/bwamem_HG004.sam ${HOME}/run_bwamem_mapping/bwamem_HG003.sam
./real_read_deeptrio_call.sh DRAGEN HG002 HG004 HG003 ${HOME}/run_dragen_mapping/dragen_output_HG002.bam ${HOME}/run_dragen_mapping/dragen_output_HG004.bam ${HOME}/run_dragen_mapping/dragen_output_HG003.bam
./real_read_deeptrio_call.sh GIRAFFE_1000GP HG002 HG004 HG003 ${HOME}/run_giraffe_mapping/giraffe_HG002.bam ${HOME}/run_giraffe_mapping/giraffe_HG004.bam ${HOME}/run_giraffe_mapping/giraffe_HG003.bam
./real_read_deeptrio_call.sh GIRAFFE_PARENTAL HG002 HG004 HG003 ${HOME}/run_giraffe_pedigree_mapping/HG002_default_vg_pedigree_outstore/HG002_merged.bam ${HOME}/run_giraffe_pedigree_mapping/HG002_default_vg_pedigree_outstore/HG004_merged.bam ${HOME}/run_giraffe_pedigree_mapping/HG002_default_vg_pedigree_outstore/HG003_merged.bam

./real_read_deeptrio_call.sh BWAMEM HG005 HG007 HG006 ${HOME}/run_bwamem_mapping/bwamem_HG005.sam ${HOME}/run_bwamem_mapping/bwamem_HG007.sam ${HOME}/run_bwamem_mapping/bwamem_HG006.sam
./real_read_deeptrio_call.sh DRAGEN HG005 HG007 HG006 ${HOME}/run_dragen_mapping/dragen_output_HG005.bam ${HOME}/run_dragen_mapping/dragen_output_HG007.bam ${HOME}/run_dragen_mapping/dragen_output_HG006.bam
./real_read_deeptrio_call.sh GIRAFFE_1000GP HG005 HG007 HG006 ${HOME}/run_giraffe_mapping/giraffe_HG005.bam ${HOME}/run_giraffe_mapping/giraffe_HG007.bam ${HOME}/run_giraffe_mapping/giraffe_HG006.bam
./real_read_deeptrio_call.sh GIRAFFE_PARENTAL HG005 HG007 HG006 ${HOME}/run_giraffe_pedigree_mapping/HG005_default_vg_pedigree_outstore/HG005_merged.bam ${HOME}/run_giraffe_pedigree_mapping/HG005_default_vg_pedigree_outstore/HG007_merged.bam ${HOME}/run_giraffe_pedigree_mapping/HG005_default_vg_pedigree_outstore/HG006_merged.bam
```

Following that, run the calling on the trained deeptrio experimental BAMs of the 1000GP-based and Parental-based Giraffe alignments.

```
./real_read_deeptrio_call.sh GIRAFFE_1000GP_TRAINED_CHR20 HG002 HG004 HG003 ${HOME}/run_giraffe_mapping/giraffe_HG002.bam ${HOME}/run_giraffe_mapping/giraffe_HG004.bam ${HOME}/run_giraffe_mapping/giraffe_HG003.bam false false true
./real_read_deeptrio_call.sh GIRAFFE_PARENTAL_TRAINED_CHR20 HG002 HG004 HG003 ${HOME}/run_giraffe_pedigree_mapping/HG002_default_vg_pedigree_outstore/HG002_merged.bam ${HOME}/run_giraffe_pedigree_mapping/HG002_default_vg_pedigree_outstore/HG004_merged.bam ${HOME}/run_giraffe_pedigree_mapping/HG002_default_vg_pedigree_outstore/HG003_merged.bam false false true

./real_read_deeptrio_call.sh GIRAFFE_1000GP_TRAINED_CHR20 HG005 HG007 HG006 ${HOME}/run_giraffe_mapping/giraffe_HG005.bam ${HOME}/run_giraffe_mapping/giraffe_HG007.bam ${HOME}/run_giraffe_mapping/giraffe_HG006.bam false false true
./real_read_deeptrio_call.sh GIRAFFE_PARENTAL_TRAINED_CHR20 HG005 HG007 HG006 ${HOME}/run_giraffe_pedigree_mapping/HG005_default_vg_pedigree_outstore/HG005_merged.bam ${HOME}/run_giraffe_pedigree_mapping/HG005_default_vg_pedigree_outstore/HG007_merged.bam ${HOME}/run_giraffe_pedigree_mapping/HG005_default_vg_pedigree_outstore/HG006_merged.bam false false true
```

Finally run the variant call evaluations on the called VCFs
```
./real_read_variant_call_evaluation.sh HG002 ${HOME}/run_dragen_genotyping/bwamem_HG002.sam_dragen_run.vcf.gz BWAMEM DRAGEN
./real_read_variant_call_evaluation.sh HG002 ${HOME}/run_dragen_genotyping/dragen_output_HG002.bam_dragen_run.vcf.gz DRAGEN DRAGEN
./real_read_variant_call_evaluation.sh HG002 ${HOME}/run_dragen_genotyping/giraffe_HG002.bam_dragen_run.vcf.gz VG_1000GP DRAGEN
./real_read_variant_call_evaluation.sh HG002 ${HOME}/run_dragen_genotyping/HG002_merged.bam_dragen_run.vcf.gz VG_PARENTAL DRAGEN

./real_read_variant_call_evaluation.sh HG005 ${HOME}/run_dragen_genotyping/bwamem_HG005.sam_dragen_run.vcf.gz BWAMEM DRAGEN
./real_read_variant_call_evaluation.sh HG005 ${HOME}/run_dragen_genotyping/dragen_output_HG005.bam_dragen_run.vcf.gz DRAGEN DRAGEN
./real_read_variant_call_evaluation.sh HG005 ${HOME}/run_dragen_genotyping/giraffe_HG005.bam_dragen_run.vcf.gz VG_1000GP DRAGEN
./real_read_variant_call_evaluation.sh HG005 ${HOME}/run_dragen_genotyping/HG005_merged.bam_dragen_run.vcf.gz VG_PARENTAL DRAGEN

./real_read_variant_call_evaluation.sh HG002 ${HOME}/run_deeptrio_genotyping/HG002_DEEPTRIO.abra_gatk_targets.BWAMEM.vcf.gz BWAMEM DEEPTRIO_DEFAULT true
./real_read_variant_call_evaluation.sh HG002 ${HOME}/run_deeptrio_genotyping/HG002_DEEPTRIO.abra_gatk_targets.DRAGEN.vcf.gz DRAGEN DEEPTRIO_DEFAULT true
./real_read_variant_call_evaluation.sh HG002 ${HOME}/run_deeptrio_genotyping/HG002_DEEPTRIO.abra_gatk_targets.GIRAFFE_1000GP.vcf.gz VG_1000GP DEEPTRIO_DEFAULT true
./real_read_variant_call_evaluation.sh HG002 ${HOME}/run_deeptrio_genotyping/HG002_DEEPTRIO.abra_gatk_targets.GIRAFFE_PARENTAL.vcf.gz VG_PARENTAL DEEPTRIO_DEFAULT true

./real_read_variant_call_evaluation.sh HG005 ${HOME}/run_deeptrio_genotyping/HG005_DEEPTRIO.abra_gatk_targets.BWAMEM.vcf.gz BWAMEM DEEPTRIO_DEFAULT true
./real_read_variant_call_evaluation.sh HG005 ${HOME}/run_deeptrio_genotyping/HG005_DEEPTRIO.abra_gatk_targets.DRAGEN.vcf.gz DRAGEN DEEPTRIO_DEFAULT true
./real_read_variant_call_evaluation.sh HG005 ${HOME}/run_deeptrio_genotyping/HG005_DEEPTRIO.abra_gatk_targets.GIRAFFE_1000GP.vcf.gz VG_1000GP DEEPTRIO_DEFAULT true
./real_read_variant_call_evaluation.sh HG005 ${HOME}/run_deeptrio_genotyping/HG005_DEEPTRIO.abra_gatk_targets.GIRAFFE_PARENTAL.vcf.gz VG_PARENTAL DEEPTRIO_DEFAULT true

./real_read_variant_call_evaluation.sh HG002 ${HOME}/run_deeptrio_genotyping/HG002_DEEPTRIO.abra_gatk_targets.BWAMEM.vcf.gz BWAMEM DEEPTRIO_DEFAULT true true
./real_read_variant_call_evaluation.sh HG002 ${HOME}/run_deeptrio_genotyping/HG002_DEEPTRIO.abra_gatk_targets.DRAGEN.vcf.gz DRAGEN DEEPTRIO_DEFAULT true true
./real_read_variant_call_evaluation.sh HG002 ${HOME}/run_deeptrio_genotyping/HG002_DEEPTRIO.abra_gatk_targets.GIRAFFE_1000GP_TRAINED_CHR20.vcf.gz VG_1000GP DEEPTRIO_TRAINED true true
./real_read_variant_call_evaluation.sh HG002 ${HOME}/run_deeptrio_genotyping/HG002_DEEPTRIO.abra_gatk_targets.GIRAFFE_PARENTAL_TRAINED_CHR20.vcf.gz VG_PARENTAL DEEPTRIO_TRAINED true true
./real_read_variant_call_evaluation.sh HG002 ${HOME}/run_deeptrio_genotyping/HG002_DEEPTRIO.abra_gatk_targets.GIRAFFE_1000GP.vcf.gz VG_1000GP DEEPTRIO_DEFAULT true true
./real_read_variant_call_evaluation.sh HG002 ${HOME}/run_deeptrio_genotyping/HG002_DEEPTRIO.abra_gatk_targets.GIRAFFE_PARENTAL.vcf.gz VG_PARENTAL DEEPTRIO_DEFAULT true true

./real_read_variant_call_evaluation.sh HG005 ${HOME}/run_deeptrio_genotyping/HG005_DEEPTRIO.abra_gatk_targets.BWAMEM.vcf.gz BWAMEM DEEPTRIO_DEFAULT true true
./real_read_variant_call_evaluation.sh HG005 ${HOME}/run_deeptrio_genotyping/HG005_DEEPTRIO.abra_gatk_targets.DRAGEN.vcf.gz DRAGEN DEEPTRIO_DEFAULT true true
./real_read_variant_call_evaluation.sh HG005 ${HOME}/run_deeptrio_genotyping/HG005_DEEPTRIO.abra_gatk_targets.GIRAFFE_1000GP_TRAINED_CHR20.vcf.gz VG_1000GP DEEPTRIO_TRAINED true true
./real_read_variant_call_evaluation.sh HG005 ${HOME}/run_deeptrio_genotyping/HG005_DEEPTRIO.abra_gatk_targets.GIRAFFE_PARENTAL_TRAINED_CHR20.vcf.gz VG_PARENTAL DEEPTRIO_TRAINED true true
./real_read_variant_call_evaluation.sh HG005 ${HOME}/run_deeptrio_genotyping/HG005_DEEPTRIO.abra_gatk_targets.GIRAFFE_1000GP.vcf.gz VG_1000GP DEEPTRIO_DEFAULT true true
./real_read_variant_call_evaluation.sh HG005 ${HOME}/run_deeptrio_genotyping/HG005_DEEPTRIO.abra_gatk_targets.GIRAFFE_PARENTAL.vcf.gz VG_PARENTAL DEEPTRIO_DEFAULT true true
```

Render the calculated roc curves into plots
```
./real_read_variant_call_evaluation.ROC_plots.sh
```

## Output files

All outputs will be located in the following directory `${HOME}/run_genotyping`.
DRAGEN-based caller outputs will be located in files named `${ALIGNER}_${SAMPLE}.*am_dragen_run.vcf.gz`
where `${ALIGNER}` is one of `bwamem`, `dragen`, `giraffe`, or `${SAMPLE}_merged`
and `${SAMPLE}` is one of `HG002` or `HG005`.


DeepTrio-based caller outputs will be located in files named `${SAMPLE}_DEEPTRIO.abra_gatk_targets.${ALIGNER}.vcf.gz`
where `${ALIGNER}` is one of `BWAMEM`, `DRAGEN`, `GIRAFFE_1000GP`, `GIRAFFE_1000GP_TRAINED_CHR20`, `GIRAFFE_PARENTAL`, or `GIRAFFE_PARENTAL_TRAINED_CHR20`
and `${SAMPLE}` is one of `HG002` or `HG005`.

Variant calling evaluation outputs will be located in files named `${EVAL}_vcfeval_output_${SAMPLE}_${REGION}_${ALIGNER}_${CALLER}${CHR20_FLAG}`
where `${EVAL}` is one of `happy` or `rtg`
and `${SAMPLE}` is one of `HG002` or `HG005`
and `${REGION}` is one of `allhighconfregions`, `alldifficultregions`, `alllowmapandsegdupregions`, `MHC`, or `CMRG`
and `${ALIGNER}` is one of `BWAMEM`, `DRAGEN`, `VG_1000GP`, or `VG_PARENTAL`
and `${CALLER}` is one of `DRAGEN`, `DEEPTRIO_DEFAULT`, `DEEPTRIO_TRAINED_chr20`.




