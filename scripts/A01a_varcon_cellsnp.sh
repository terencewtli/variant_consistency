#!/bin/bash
#$ -cwd
#$ -o logs/A01a_pileup.$JOB_ID.$TASK_ID
#$ -j y
#$ -N A01a_pileup
#$ -l h_data=1G,h_rt=1:00:00
#$ -pe shared 1
#$ -t 1-12:1

echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "

source ~/.bashrc

PROJDIR=/u/home/t/terencew/project-cluo/igvf/50_line/multiome/rna
cd $PROJDIR

TEMPLATE=/u/project/cluo/terencew/demux_benchmark/template_demux/

CR_BASE=$PROJDIR/mapping/cr_arc/

# ID=1
ID=${SGE_TASK_ID}
SAMPLES=$PROJDIR/txt/mapping/batch2/gex_samples.txt
SAMPLE=$(head -${ID} $SAMPLES | tail -1)

INDIR=$PROJDIR/demux/wgs/
cd $INDIR

GEX_BAM=${CR_BASE}/$SAMPLE/outs/gex_possorted_bam.bam
ATAC_BAM=${CR_BASE}/$SAMPLE/outs/atac_possorted_bam.bam
BARCODES=${CR_BASE}/$SAMPLE/outs/filtered_feature_bc_matrix/barcodes.tsv.gz

VCF=$PROJDIR/vcf/wgs/snp_filt2_demux.rm_missing.vcf.gz
VARCON_GEX=$INDIR/varcon_cellsnp/gex/$SAMPLE/
VARCON_ATAC=$INDIR/varcon_cellsnp/atac/$SAMPLE/
qsub $TEMPLATE/varcon_cellsnp_gex.sh $GEX_BAM $BARCODES $VCF $VARCON_GEX
qsub $TEMPLATE/varcon_cellsnp_atac.sh $ATAC_BAM $BARCODES $VCF $VARCON_ATAC

echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "
