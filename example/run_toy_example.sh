#!/bin/bash
### Variant consistency tutorial on a synthetic toy experiment (runs in seconds; needs numpy, pandas, scipy).
### run from the repository root: bash example/run_toy_example.sh
set -e
OUT=example/toy

### 1. make the toy data: cellSNP-lite style pileup, donor genotypes, donor list, demultiplexing assignments
python example/make_toy_example.py --outdir $OUT

### 2. variant consistency counts per droplet and donor (C1, C2, I1, I2) at SNP depth thresholds 0 and 10
python varcon_fast.py -c $OUT/cellsnp -i $OUT/genotypes.vcf.gz -d $OUT/donors.txt -o $OUT/varcon --covs 0,10

### 3. per-method tables: counts at each droplet's assigned donor, for every method column in demux.tsv
mkdir -p $OUT/varcon/per_method
python 04_organize_varcon.py -i $OUT/varcon -d $OUT/donors.txt -c 0 -v $OUT/demux.tsv -o $OUT/varcon/per_method

### 4. misassigned droplets carry more I1 (alleles that match another donor but not the assigned one)
python example/summarize_toy.py $OUT/varcon/per_method/toy_method_0.csv $OUT/toy_method_status.tsv
