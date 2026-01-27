#!/bin/bash
#$ -cwd
#$ -o logs/A01d_get_varcon_counts.$JOB_ID.$TASK_ID
#$ -j y
#$ -N A01d_get_varcon_counts
#$ -l h_data=2G,h_rt=2:00:00
#$ -pe shared 4
#$ -t 1-3:1

echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "

source ~/.bashrc

conda activate scFates

PROJDIR=/u/home/t/terencew/project-cluo/igvf/pilot/multiome
cd $PROJDIR

TEMPLATE=/u/home/t/terencew/project-cluo/demux_benchmark/template_demux/

# ID=1
ID=${SGE_TASK_ID}
SAMPLES=$PROJDIR/txt/samples.txt
SAMPLE=$(head -${ID} $SAMPLES | tail -1)

N_DONORS=4

INDIR=$PROJDIR/demux/regular

GEX_VCF=$PROJDIR/vcf/pilot_000.vcf.gz
ATAC_VCF=$PROJDIR/vcf/pilot_000.vcf.gz
DONORS=$PROJDIR/txt/donors.txt

VARCON=/u/home/t/terencew/project-cluo/programs/var_consistency/03_count_varcon.py

# COVS=(0 5 10 15 20)
COVS=(0 10 20)

GEX_INDIR=$INDIR/var_consistency/gex/$SAMPLE/consistency/
GEX_PKLDIR=$INDIR/var_consistency/gex/$SAMPLE/dict/

time for COV in ${COVS[@]};
do
    GEX_OUTDIR=$INDIR/var_consistency/gex/$SAMPLE/csv/cov_${COV}
    python $VARCON -i $GEX_INDIR -p $GEX_PKLDIR \
        -d $DONORS -o $GEX_OUTDIR -c $COV
done

ATAC_INDIR=$INDIR/var_consistency/atac/$SAMPLE/consistency/
ATAC_PKLDIR=$INDIR/var_consistency/atac/$SAMPLE/dict/

time for COV in ${COVS[@]};
do
    ATAC_OUTDIR=$INDIR/var_consistency/atac/$SAMPLE/csv/cov_${COV}
    python $VARCON -i $ATAC_INDIR -p $ATAC_PKLDIR \
        -d $DONORS -o $ATAC_OUTDIR -c $COV
done

echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `date `
echo " "
