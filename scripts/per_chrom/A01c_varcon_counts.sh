#!/bin/bash
#$ -cwd
#$ -o logs/A01c_varcon_counts.$JOB_ID
#$ -j y
#$ -N A01c_varcon_counts
#$ -l h_data=1G,h_rt=4:00:00
#$ -pe shared 12

echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "

source ~/.bashrc

conda activate scFates

### this means that the script reads the first argument as
### chromosome, second as sample
CHROM=$1
SAMPLE=$2

DONORS=$PROJDIR/txt/donors.txt
INDIR=$PROJDIR/demux/wgs

TEMPLATE=/u/home/t/terencew/project-cluo/demux_benchmark/template_demux/
VARCON=/u/project/cluo/terencew/programs/variant_consistency/01_con_counts_multithread.py

GEX_VCF=$PROJDIR/vcf/wgs/by_chrom/snp_filt2_demux.rm_missing.${CHROM}.vcf.gz
GEX_CELLSNP=$INDIR/varcon_cellsnp/gex/$SAMPLE/by_chrom/$CHROM
GEX_OUTDIR=$INDIR/var_consistency/gex/$SAMPLE/by_chrom/$CHROM/consistency

time python $VARCON -c $GEX_CELLSNP \
    -i $GEX_VCF -d $DONORS -o $GEX_OUTDIR -t 12

echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `date `
echo " "
