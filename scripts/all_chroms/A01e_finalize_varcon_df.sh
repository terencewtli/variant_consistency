#!/bin/bash
#$ -cwd
#$ -o logs/A01e_finalize_varcon.$JOB_ID.$TASK_ID
#$ -j y
#$ -N A01e_finalize_varcon
#$ -l h_data=2G,h_rt=4:00:00
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
DONORS=$PROJDIR/txt/donors.txt

GEX_VCF=$PROJDIR/vcf/pilot_000.vcf.gz
ATAC_VCF=$PROJDIR/vcf/pilot_000.vcf.gz

VARCON=/u/home/t/terencew/project-cluo/programs/var_consistency/04_organize_varcon.py

INDIR=$PROJDIR/demux/regular

# COVS=(0 5 10 15 20)
COVS=(0 10 20)

GEX_DEMUX=$PROJDIR/csv/demux/inter_gex.csv

time for COV in ${COVS[@]};
do
    GEX_INDIR=$INDIR/var_consistency/gex/$SAMPLE/csv/cov_${COV}
    GEX_OUTDIR=$INDIR/var_consistency/gex/$SAMPLE/csv/cov_${COV}
    python $VARCON -i $GEX_INDIR -d $DONORS \
        -v $GEX_DEMUX -o $GEX_OUTDIR -c $COV
done

ATAC_DEMUX=$PROJDIR/csv/demux/inter_atac.csv

time for COV in ${COVS[@]};
do
    ATAC_INDIR=$INDIR/var_consistency/atac/$SAMPLE/csv/cov_${COV}
    ATAC_OUTDIR=$INDIR/var_consistency/atac/$SAMPLE/csv/cov_${COV}
    python $VARCON -i $ATAC_INDIR -d $DONORS \
        -v $ATAC_DEMUX -o $ATAC_OUTDIR -c $COV
done

echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `date `
echo " "
