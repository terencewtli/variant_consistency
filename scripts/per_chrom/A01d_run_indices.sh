#!/bin/bash
#$ -cwd
#$ -o logs/A01d_run_indices.$JOB_ID.$TASK_ID
#$ -j y
#$ -N A01d_run_indices
#$ -l h_data=1G,h_rt=1:00:00
#$ -pe shared 1
#$ -t 1-51:1

echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "

source ~/.bashrc

conda activate scFates

PROJDIR=/u/home/t/terencew/project-cluo/igvf/50_line/multiome/rna # <-
cd $PROJDIR

# ID=1
ID=${SGE_TASK_ID}
SAMPLES=$PROJDIR/txt/mapping/batch2/gex_samples.txt
SAMPLE=$(head -${ID} $SAMPLES | tail -1)

VARCON=/u/project/cluo/terencew/demux_benchmark/programs/for_nadia/50line_chrom/per_chrom/A01d_varcon_indices.sh

for c in `seq 1 22`;
do
    CHROM=chr$c
    qsub $VARCON $CHROM $SAMPLE
done

echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "
