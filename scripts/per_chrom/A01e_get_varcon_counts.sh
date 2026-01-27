#!/bin/bash
#$ -cwd
#$ -o logs/A01e_get_varcon_counts.$JOB_ID.$TASK_ID
#$ -j y
#$ -N A01e_get_varcon_counts
#$ -l h_data=1G,h_rt=4:00:00
#$ -pe shared 4
#$ -t 1-51:1

echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "

source ~/.bashrc

conda activate scFates

INDIR=$PROJDIR/demux/wgs
cd $INDIR

CHROM=$1
SAMPLE=$2

VARCON=/u/home/t/terencew/project-cluo/programs/variant_consistnecy/03_count_varcon.py

COV=0

GEX_VCF=$PROJDIR/vcf/wgs/by_chrom/snp_filt2_demux.rm_missing.${CHROM}.vcf.gz

GEX_INDIR=$PROJDIR/demux/wgs/var_consistency/gex/$SAMPLE/by_chrom/$CHROM/consistency/ # <-
GEX_PKLDIR=$PROJDIR/demux/wgs/var_consistency/gex/$SAMPLE/by_chrom/$CHROM/dict/
GEX_OUTDIR=$PROJDIR/demux/wgs/var_consistency/gex/$SAMPLE/by_chrom/$CHROM/csv/cov_${COV}
mkdir -p $GEX_OUTDIR

python $VARCON -i $GEX_INDIR -p $GEX_PKLDIR \
      -d $DONORS -o $GEX_OUTDIR -c $COV

echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "
