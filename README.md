# Purpose of repository

This repository contains code used to measure ambient contamination of
genotype-based demultiplexing methods using reference genotypes in
multiplexed single-cell experiments from the following manuscript:

Li, Terence, Marcus Alvarez, Cuining Liu, Kevin Abuhanna, Yu Sun, Jason Ernst, Kathrin Plath, Brunilda Balliu, Chongyuan Luo, and Noah Zaitlen. "The impact of ambient contamination on demultiplexing methods for single-nucleus multiome experiments." bioRxiv (2025).

Elife reviewed preprint:
[https://elifesciences.org/reviewed-preprints/106769
](https://elifesciences.org/reviewed-preprints/106769)

bioRxiv:
[https://www.biorxiv.org/content/10.1101/2025.02.06.636969v1.full.pdf
](https://www.biorxiv.org/content/10.1101/2025.02.06.636969v1.full.pdf)


# What variant consistency measures

For every droplet that a demultiplexing method assigns to a donor, variant consistency compares the alleles observed in
that droplet with the reference genotypes of all multiplexed donors. Each allele count falls into one of four mutually
exclusive categories, relative to the assigned donor *d*:

| category | allele is consistent with | interpretation |
|---|---|---|
| **C1** | *d* only | informative support for the assignment |
| **C2** | *d* and at least one other donor | uninformative (shared alleles) |
| **I1** | another donor but not *d* | reads from another donor: ambient contamination, a missed doublet or a wrong assignment (plus a small error background) |
| **I2** | no donor | sequencing, alignment or genotyping error |

Rates are each category divided by C1 + C2 + I1 + I2. Comparing I1 rates between methods on the same data shows which
method retains more contaminated or misassigned droplets.

# Installation

The counting step needs Python 3 with numpy, pandas and scipy:

```
conda env create -f workflow/envs/var_con.yaml
conda activate var_con_env
```

Allele counts per droplet come from [cellsnp-lite](https://github.com/single-cell-genetics/cellsnp-lite)
(`workflow/envs/cellsnp.yaml`). The Snakemake workflow (below) creates both environments for you.

# Tutorial: a toy example

`example/` contains a fully synthetic multiplexed experiment (4 donors, 2,000 SNPs, 400 droplets with 10% doublets and
~15% ambient contamination) and a demultiplexing "method" that misassigns 20 singlets and misses 15 doublets. It runs in
seconds:

```
bash example/run_toy_example.sh
```

The script runs three steps that you would also run on real data:

```
# 1. make the toy data (cellsnp-lite style pileup, genotypes, donor list, assignments)
python example/make_toy_example.py --outdir example/toy

# 2. counts per droplet and donor for each category, at SNP depth thresholds 0 and 10
python varcon_fast.py -c example/toy/cellsnp -i example/toy/genotypes.vcf.gz \
    -d example/toy/donors.txt -o example/toy/varcon --covs 0,10

# 3. per-method tables: counts at each droplet's assigned donor, for every method column in demux.tsv
python 04_organize_varcon.py -i example/toy/varcon -d example/toy/donors.txt -c 0 \
    -v example/toy/demux.tsv -o example/toy/varcon/per_method
```

and prints the median composition of the toy method's singlet calls:

```
             droplets   C1    C2    I1   I2
status
correct           340  1.5  95.9   2.3  0.0
doublet            15  0.8  89.3  10.0  0.0
wrong_donor        20  0.0  83.3  16.5  0.0
```

Correctly assigned droplets carry a little I1 (from ambient contamination); missed doublets and wrong-donor calls carry
much more.

# Running on your data

## Inputs

- **Pileup:** cellsnp-lite output for one library at the genotyped sites, e.g.
  ```
  cellsnp-lite -s possorted_bam.bam -b barcodes.tsv.gz -R genotypes.vcf.gz -O cellsnp_out -p 8 \
      --minMAF 0.1 --minCOUNT 1 --gzip            # add --UMItag None for snATAC-seq
  ```
- **Genotypes:** a VCF with a GT field for every multiplexed donor (we used imputed genotypes with R² > 0.90).
- **Donor list:** a text file with one donor ID per line, matching the VCF header.
- **Assignments (for per-method tables):** a tab-separated table, one row per barcode and one column per method; values are
  a donor ID, `doublet` or `unassigned`.

## One library, by hand

```
python varcon_fast.py -c cellsnp_out -i genotypes.vcf.gz -d donors.txt -o varcon_out --covs 0,10,20
python 04_organize_varcon.py -i varcon_out -d donors.txt -c 0 -v assignments.tsv -o varcon_out/per_method
```

`varcon_fast.py` writes `{c1,c2,i1,i2}_df_{cov}.csv` (droplets × donors) for each SNP depth threshold (SNPs are kept when
their experiment-wide depth exceeds the threshold). `04_organize_varcon.py` writes `<method>_<cov>.csv` (C1, C2, I1, I2
and the assigned donor per droplet) and `inter_<cov>_varcon.csv` (pooled rates per method).

## Many libraries, with Snakemake

`workflow/Snakefile` runs cellsnp-lite, the counts and (optionally) the per-method tables for every library listed in
`config/config.yaml`:

```
conda env create -f workflow/envs/snakemake.yaml
conda activate snakemake_env
# edit config/config.yaml: genotype VCF, donor list, and one entry per library (BAM, barcodes, gex/atac, optional assignments)
snakemake -np                              # dry run
snakemake --use-conda --cores 8            # local
qsub submit_snakemake.sh                   # SGE cluster (Hoffman2 example)
```

Outputs, per library: `00_cellsnp/` (pileup), `02_varcon/` (counts), `03_per_method/` (per-method tables) and `logs/`.

# Original implementation

`01_con_counts_multithread.py`, `02_get_con_indices.py` and `03_count_varcon*.py` are the original step-by-step
implementation. `varcon_fast.py` produces identical counts in a single vectorised pass (a 10x Multiome library takes a few
minutes, mostly reading the VCF) and is what the tutorial and workflow use.

# Acknowledgements

The Snakemake workflow was written by **Nadia Shaw** ([github.com/nadiashaw0405/my_var_con](https://github.com/nadiashaw0405/my_var_con))
and adapted here to use `varcon_fast.py`, the cellsnp-lite settings of the manuscript, and optional per-method tables.
