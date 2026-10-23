#!/usr/bin/env python
import os
import sys
import gzip
import time
import argparse
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix

cats = ['c1', 'c2', 'i1', 'i2']

def read_mtx(path):
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as f:
        n_header = 0
        for line in f:
            n_header += 1
            if not line.startswith('%'):
                n_row, n_col, _ = map(int, line.split())
                break
    df = pd.read_csv(path, sep=r'\s+', skiprows=n_header, header=None, names=['r', 'c', 'v'],
                     dtype={'r': np.int64, 'c': np.int64, 'v': np.int64})
    return coo_matrix((df['v'].values, (df['r'].values - 1, df['c'].values - 1)), shape=(n_row, n_col))


def find(cellsnp_dir, name):
    for f in [name, f'{name}.gz']:
        if os.path.exists(f'{cellsnp_dir}/{f}'):
            return f'{cellsnp_dir}/{f}'
    raise FileNotFoundError(f'{cellsnp_dir}/{name}')


def vcf_header(path):
    opener = gzip.open if path.endswith('.gz') else open
    last = None
    with opener(path, 'rt') as f:
        for line in f:
            if not line.startswith('#'):
                break
            last = line
    return last.lstrip('#').rstrip('\n').split('\t')

def genotype_codes(gts):
    ### 0 = het, 1 = hom ref, 2 = hom alt, 3 = other (missing, ...): same string tests as 01's get_con_counts
    gts = pd.Series(gts).astype(str).str.split(':').str[0]
    has0 = gts.str.contains('0', regex=False).values
    has1 = gts.str.contains('1', regex=False).values
    code = np.full(len(gts), 3, dtype=np.int8)
    code[has0 & has1] = 0
    code[has0 & ~has1] = 1
    code[~has0 & has1] = 2
    return code

def load_genotypes(vcf, donors):
    ### genotype codes (sites x donors) indexed by chrom_pos; first copy of duplicated sites kept
    cols = vcf_header(vcf)
    missing = [d for d in donors if d not in cols]
    if missing:
        raise ValueError(f'donors not in VCF header: {missing}')
    geno = pd.read_csv(vcf, sep='\t', header=None, comment='#', usecols=[0, 1] + [cols.index(d) for d in donors], dtype=str)
    geno.columns = ['CHROM', 'POS'] + donors
    geno.index = geno['CHROM'] + '_' + geno['POS']
    geno = geno[~geno.index.duplicated()]
    return pd.DataFrame(np.stack([genotype_codes(geno[d].values) for d in donors], axis=1), index=geno.index, columns=donors)


def count(cellsnp_dir, geno, covs=(0, 10, 20), log=print):
    ### returns {(cat, cov): DataFrame cells x donors}
    t0 = time.time()
    donors = list(geno.columns)
    n_donors = len(donors)
    barcodes = pd.read_csv(find(cellsnp_dir, 'cellSNP.samples.tsv'), sep='\t', header=None)[0].values
    snps = pd.read_csv(find(cellsnp_dir, 'cellSNP.base.vcf'), sep='\t', header=None, comment='#', usecols=[0, 1, 7], dtype={0: str})
    snps.index = snps[0] + '_' + snps[1].astype(str)
    snps['DP'] = snps[7].str.split(';').str[1].str.split('=').str[1].astype(int)

    keep = snps.index.isin(geno.index) & ~snps.index.duplicated()
    snp_rows = np.flatnonzero(keep)
    codes = geno.loc[snps.index[keep]].values                      # n_kept x n_donors

    dp = read_mtx(find(cellsnp_dir, 'cellSNP.tag.DP.mtx')).tocsr()
    ad = read_mtx(find(cellsnp_dir, 'cellSNP.tag.AD.mtx')).tocsr()
    row_map = np.full(dp.shape[0], -1, dtype=np.int64)
    row_map[snp_rows] = np.arange(len(snp_rows))
    dp = dp.tocoo()
    m = (row_map[dp.row] >= 0) & (dp.data > 0)
    r, c, d = dp.row[m], dp.col[m], dp.data[m].astype(np.int32)
    a = np.asarray(ad[r, c]).ravel().astype(np.int32)
    k = row_map[r]
    del dp, ad
    log(f'{len(snps)} sites, {int(keep.sum())} in VCF, {len(barcodes)} cells, {n_donors} donors, {len(d)} entries '
        f'({time.time() - t0:.0f}s)')

    def consistent(j):
        g = codes[k, j]
        return np.where(g == 0, d, np.where(g == 1, d - a, np.where(g == 2, a, 0))).astype(np.int32)

    ### pass 1: how many donors each entry is consistent / inconsistent with
    n_con = np.zeros(len(d), dtype=np.int16)
    n_incon = np.zeros(len(d), dtype=np.int16)
    for j in range(n_donors):
        con = consistent(j)
        n_con += con > 0
        n_incon += (d - con) > 0
    c1_m, c2_m = n_con == 1, n_con > 1
    i1_m, i2_m = (n_incon > 0) & (n_incon < n_donors), n_incon == n_donors

    ### pass 2: per-donor category sums per cell, for each SNP depth threshold
    site_dp = snps['DP'].values[snp_rows][k]
    sels = {cov: site_dp > cov for cov in covs}
    out = {(cat, cov): np.zeros((len(barcodes), n_donors), dtype=np.int64) for cat in cats for cov in covs}
    for j in range(n_donors):
        con = consistent(j)
        incon = d - con
        w = {'c1': np.where((con > 0) & c1_m, con, 0), 'c2': np.where((con > 0) & c2_m, con, 0),
             'i1': np.where((incon > 0) & i1_m, incon, 0), 'i2': np.where((con == 0) & i2_m, incon, 0)}
        for cov, sel in sels.items():
            for cat in cats:
                out[(cat, cov)][:, j] = np.bincount(c[sel], weights=w[cat][sel], minlength=len(barcodes)).astype(np.int64)
    log(f'counted ({time.time() - t0:.0f}s)')
    return {key: pd.DataFrame(v, index=barcodes, columns=donors) for key, v in out.items()}


def main():
    ap = argparse.ArgumentParser(description='Fast variant consistency counts per cell and donor')
    ap.add_argument('-c', '--cellsnp_dir', required=True)
    ap.add_argument('-i', '--vcf', required=True, help='VCF with reference genotypes')
    ap.add_argument('-d', '--donors', required=True, help='txt file, one multiplexed donor per line')
    ap.add_argument('-o', '--outdir', required=True)
    ap.add_argument('--covs', default='0,10,20', help='experiment-wide SNP depth thresholds (kept if DP > cov)')
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    donors = list(np.loadtxt(args.donors, dtype=str, ndmin=1))
    covs = [int(x) for x in args.covs.split(',')]
    res = count(args.cellsnp_dir, load_genotypes(args.vcf, donors), covs, log=lambda x: print(x, flush=True))
    for (cat, cov), df in res.items():
        df.to_csv(f'{args.outdir}/{cat}_df_{cov}.csv', sep='\t')

if __name__ == '__main__':
    main()
