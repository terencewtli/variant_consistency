#!/bin/bash
#$ -cwd
#$ -o logs/A01c_varcon_indices.$JOB_ID.$TASK_ID
#$ -j y
#$ -N A01c_varcon_indices
#$ -l h_data=2G,h_rt=12:00:00
#$ -pe shared 8
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

INDIR=$PROJDIR/demux/regular
cd $INDIR

DONORS=$PROJDIR/txt/donors.txt

VARCON=/u/home/t/terencew/project-cluo/programs/var_consistency/02_get_con_indices.py

GEX_INDIR=$INDIR/var_consistency/gex/$SAMPLE/consistency/
GEX_OUTDIR=$INDIR/var_consistency/gex/$SAMPLE/dict/

time python $VARCON -i $GEX_INDIR -d $DONORS \
      -o $GEX_OUTDIR -t 8

ATAC_INDIR=$INDIR/var_consistency/atac/$SAMPLE/consistency/
ATAC_OUTDIR=$INDIR/var_consistency/atac/$SAMPLE/dict/

time python $VARCON -i $ATAC_INDIR -d $DONORS \
      -o $ATAC_OUTDIR -t 8

echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `date `
echo " "
