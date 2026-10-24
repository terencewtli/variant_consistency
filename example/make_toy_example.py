#!/usr/bin/env python
import os
import gzip
import argparse
import numpy as np
import pandas as pd

### Toy multiplexed experiment for the variant consistency tutorial (fully synthetic: no real genotypes or reads).
### 4 donors, 2,000 biallelic SNPs on one chromosome, 400 droplets (360 singlets + 40 doublets). Each droplet covers ~150 SNPs
### with 1-3 reads; a fraction of reads comes from the ambient pool (all donors), plus a 1% sequencing error.
### Writes, in --outdir:
###   cellsnp/              cellSNP-lite style output (cellSNP.base.vcf.gz, cellSNP.samples.tsv, cellSNP.tag.{AD,DP,OTH}.mtx)
###   genotypes.vcf.gz      donor genotypes (GT)
###   donors.txt            donor IDs, one per line
###   demux.tsv             droplet assignments: "truth" and "toy_method" (truth with 20 singlets swapped to another donor
###                         and 15 doublets called singlets), in the format 04_organize_varcon.py expects
### usage: python make_toy_example.py --outdir toy

rng = np.random.default_rng(1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--outdir', default='toy')
    ap.add_argument('--n_snps', type=int, default=2000)
    ap.add_argument('--n_singlets', type=int, default=360)
    ap.add_argument('--n_doublets', type=int, default=40)
    ap.add_argument('--ambient', type=float, default=0.15, help='mean fraction of reads from the ambient pool')
    args = ap.parse_args()
    out = args.outdir
    os.makedirs(f'{out}/cellsnp', exist_ok=True)

    donors = ['donorA', 'donorB', 'donorC', 'donorD']
    n_donors = len(donors)
    pos = np.sort(rng.choice(np.arange(1_000_000, 50_000_000), args.n_snps, replace=False))
    maf = rng.uniform(0.1, 0.5, args.n_snps)
    geno = rng.binomial(2, maf[:, None], size=(args.n_snps, n_donors))     # alt allele dosage 0/1/2 per donor
    alt_frac = geno / 2                                                   # expected alt fraction of a donor's reads

    ### droplets: singlets of one donor, doublets of two different donors
    n_cells = args.n_singlets + args.n_doublets
    barcodes = [f'TOY{i:04d}-1' for i in range(n_cells)]
    first = rng.integers(0, n_donors, n_cells)
    second = (first + rng.integers(1, n_donors, n_cells)) % n_donors
    is_doublet = np.arange(n_cells) >= args.n_singlets
    ambient = np.clip(rng.normal(args.ambient, 0.05, n_cells), 0, 0.6)
    pool = alt_frac.mean(axis=1)                                          # ambient = equal mix of all donors

    rows, cols, ad, dp = [], [], [], []
    for cidx in range(n_cells):
        snps = rng.choice(args.n_snps, rng.integers(100, 200), replace=False)
        depth = rng.integers(1, 4, len(snps))
        native = alt_frac[snps, first[cidx]]
        if is_doublet[cidx]:
            native = (native + alt_frac[snps, second[cidx]]) / 2
        p_alt = (1 - ambient[cidx]) * native + ambient[cidx] * pool[snps]
        p_alt = p_alt * 0.99 + (1 - p_alt) * 0.01                         # 1% sequencing error
        alt = rng.binomial(depth, p_alt)
        rows += list(snps)
        cols += [cidx] * len(snps)
        ad += list(alt)
        dp += list(depth)
    mtx = pd.DataFrame({'snp': rows, 'cell': cols, 'ad': ad, 'dp': dp})

    ### cellSNP-lite style output (SNP x cell MatrixMarket, 1-based; INFO = AD;DP;OTH)
    site = mtx.groupby('snp')[['ad', 'dp']].sum().reindex(range(args.n_snps), fill_value=0)
    with gzip.open(f'{out}/cellsnp/cellSNP.base.vcf.gz', 'wt') as f:
        f.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for i in range(args.n_snps):
            f.write(f'chr1\t{pos[i]}\t.\tA\tG\t.\tPASS\tAD={site.ad[i]};DP={site.dp[i]};OTH=0\n')
    pd.Series(barcodes).to_csv(f'{out}/cellsnp/cellSNP.samples.tsv', index=False, header=False)
    for tag, col in [('AD', 'ad'), ('DP', 'dp'), ('OTH', None)]:
        x = mtx[mtx[col] > 0] if col else mtx.iloc[:0]
        with open(f'{out}/cellsnp/cellSNP.tag.{tag}.mtx', 'w') as f:
            f.write('%%MatrixMarket matrix coordinate integer general\n%\n')
            f.write(f'{args.n_snps}\t{n_cells}\t{len(x)}\n')
            for r in x.itertuples(index=False):
                f.write(f'{r.snp + 1}\t{r.cell + 1}\t{getattr(r, col) if col else 0}\n')

    ### donor genotypes
    gt = {0: '0|0', 1: '0|1', 2: '1|1'}
    with gzip.open(f'{out}/genotypes.vcf.gz', 'wt') as f:
        f.write('##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(donors) + '\n')
        for i in range(args.n_snps):
            f.write(f'chr1\t{pos[i]}\t.\tA\tG\t.\tPASS\t.\tGT\t' + '\t'.join(gt[g] for g in geno[i]) + '\n')
    pd.Series(donors).to_csv(f'{out}/donors.txt', index=False, header=False)

    ### assignments: truth, and a toy method with wrong donors and missed doublets
    truth = np.where(is_doublet, 'doublet', np.array(donors)[first])
    method = truth.copy()
    sng = np.flatnonzero(~is_doublet)
    swap = rng.choice(sng, 20, replace=False)
    method[swap] = np.array(donors)[(first[swap] + 1) % n_donors]
    dbl = rng.choice(np.flatnonzero(is_doublet), 15, replace=False)
    method[dbl] = np.array(donors)[first[dbl]]
    pd.DataFrame({'truth': truth, 'toy_method': method}, index=barcodes).to_csv(f'{out}/demux.tsv', sep='\t')
    pd.DataFrame({'status': np.where(is_doublet, 'doublet', np.where(np.isin(np.arange(n_cells), swap), 'wrong_donor', 'correct'))},
                 index=barcodes).to_csv(f'{out}/toy_method_status.tsv', sep='\t')
    print(f'wrote {out}/: {n_cells} droplets ({args.n_doublets} doublets), {args.n_snps} SNPs, {len(mtx)} pileup entries')


if __name__ == '__main__':
    main()
