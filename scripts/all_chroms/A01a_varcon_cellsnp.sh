#!/bin/bash
#$ -cwd
#$ -o logs/A1a_varcon_cellsnp.$JOB_ID.$TASK_ID
#$ -j y
#$ -N A1a_varcon_cellsnp
#$ -l h_data=1G,h_rt=1:00:00
#$ -pe shared 1
#$ -t 1-3:1

echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "

source ~/.bashrc

conda activate scFates

PROJDIR=/u/home/t/terencew/project-cluo/igvf/pilot/multiome
cd $PROJDIR

TEMPLATE=/u/project/cluo/terencew/demux_benchmark/template_demux/

# ID=1
ID=${SGE_TASK_ID}
SAMPLES=$PROJDIR/txt/samples.txt
SAMPLE=$(head -${ID} $SAMPLES | tail -1)

N_DONORS=4
DONORS=$PROJDIR/txt/donors.txt

CR_BASE=$PROJDIR/mapping/cr_arc/igvf_ref
GEX_BAM=$CR_BASE/$SAMPLE/outs/gex_possorted_bam.bam
ATAC_BAM=$CR_BASE/$SAMPLE/outs/atac_possorted_bam.bam
BARCODES=${CR_BASE}/$SAMPLE/outs/filtered_feature_bc_matrix/barcodes.tsv.gz

VARCON=/u/home/t/terencew/project-cluo/programs/var_consistency/01_con_counts_multithread.py

DEMUX_DIR=$PROJDIR/demux/regular
GENO_VCF=$PROJDIR/vcf/pilot_000.vcf.gz
VARCON_GEX=$DEMUX_DIR/varcon_cellsnp/gex/$SAMPLE/
VARCON_ATAC=$DEMUX_DIR/varcon_cellsnp/atac/$SAMPLE/
qsub $TEMPLATE/varcon_cellsnp_gex.sh $GEX_BAM $BARCODES $GENO_VCF $VARCON_GEX
qsub $TEMPLATE/varcon_cellsnp_atac.sh $ATAC_BAM $BARCODES $GENO_VCF $VARCON_ATAC

echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "
