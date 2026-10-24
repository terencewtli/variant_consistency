#!/bin/bash
#$ -cwd
#$ -j y
#$ -o logs/snakemake_main.log
#$ -l h_rt=24:00:00,h_data=8G
#$ -N varcon_snakemake

### Run the workflow on an SGE cluster (Hoffman2 example; adapted from github.com/nadiashaw0405/my_var_con).
### Create the base environment once: conda env create -f workflow/envs/snakemake.yaml
. /u/local/Modules/default/init/modules.sh
module load anaconda3
source /u/local/apps/anaconda3/2020.11/etc/profile.d/conda.sh
conda activate snakemake_env
mkdir -p logs/cluster

snakemake --unlock
time snakemake \
    --use-conda \
    --executor cluster-generic \
    --cluster-generic-submit-cmd "qsub -l h_data=4G,h_rt=12:00:00 -pe shared {threads} -cwd -j y -o logs/cluster/ -V" \
    --jobs 4 \
    --default-resources \
    --latency-wait 60 \
    --rerun-incomplete \
    --printshellcmds \
    --cores 12
