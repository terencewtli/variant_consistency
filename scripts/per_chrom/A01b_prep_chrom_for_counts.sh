#!/bin/bash
#$ -cwd
#$ -o logs/A00_prep_cellsnp_chrom.$JOB_ID.$TASK_ID
#$ -j y
#$ -N A00_prep_cellsnp_chrom
#$ -l h_data=2G,h_rt=2:00:00
#$ -pe shared 2
#$ -t 1-51:1

echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "

source ~/.bashrc

conda activate scFates

PROJDIR=/u/home/t/terencew/project-cluo/igvf/50_line/multiome/rna
cd $PROJDIR

TEMPLATE=/u/home/t/terencew/project-cluo/demux_benchmark/template_demux

# ID=1
ID=${SGE_TASK_ID}
SAMPLES=$PROJDIR/txt/mapping/batch2/gex_samples.txt
SAMPLE=$(head -${ID} $SAMPLES | tail -1)

SUBSET=$TEMPLATE/by_chrom/subset_cellsnp.py
echo $SUBSET

### mkdir parent directory
### note: change $PROJDIR/demux/wgs/var_consistency to your parent directory
INDIR=$PROJDIR/demux/wgs/var_consistency/gex/$SAMPLE/
mkdir -p $INDIR/by_chrom

INDIR=$PROJDIR/demux/wgs/var_consistency/atac/$SAMPLE/
mkdir -p $INDIR/by_chrom

for c in `seq 1 22`;
do
    CHROM=chr${c}
    INDIR=$PROJDIR/demux/wgs/varcon_cellsnp/gex/${SAMPLE}
    python $SUBSET $INDIR $CHROM

    CHROM=chr${c}
    INDIR=$PROJDIR/demux/wgs/varcon_cellsnp/atac/${SAMPLE}
    python $SUBSET $INDIR $CHROM

    ### mkdir per chrom
    INDIR=$PROJDIR/demux/wgs/var_consistency/gex/$SAMPLE/by_chrom/
    mkdir -p $INDIR/$CHROM \
        $INDIR/$CHROM/consistency \
        $INDIR/$CHROM/csv \
        $INDIR/$CHROM/dict \

    INDIR=$PROJDIR/demux/wgs/var_consistency/gex/$SAMPLE/by_chrom/
    mkdir -p $INDIR/$CHROM \
        $INDIR/$CHROM/consistency \
        $INDIR/$CHROM/csv \
        $INDIR/$CHROM/dict \
done

echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "
