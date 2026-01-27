#!/bin/bash
#$ -cwd
#$ -o logs/A01d_varcon_indices.$JOB_ID.$TASK_ID
#$ -j y
#$ -N A01d_varcon_indices
#$ -l h_data=1G,h_rt=2:00:00
#$ -pe shared 12
#$ -t 1-51:1

echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "

source ~/.bashrc

conda activate scFates # <-

TEMPLATE=/u/home/t/terencew/project-cluo/demux_benchmark/template_demux

INDIR=$PROJDIR/demux/wgs # <-
cd $INDIR

CHROM=$1
SAMPLE=$2

DONORS=$PROJDIR/txt/donors.txt
VARCON=/u/home/t/terencew/project-cluo/programs/variant_consistency/02_get_con_indices.py

GEX_INDIR=$INDIR/var_consistency/gex/$SAMPLE/by_chrom/$CHROM/consistency/
GEX_OUTDIR=$INDIR/var_consistency/gex/$SAMPLE/by_chrom/$CHROM/dict/

python $VARCON -i $GEX_INDIR -d $DONORS \
      -o $GEX_OUTDIR -t 12

echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `date `
echo " "
