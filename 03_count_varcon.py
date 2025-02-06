#!/usr/bin/env python

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, vstack, hstack
from scipy.io import mmread, mmwrite

import argparse
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import pickle

parser = argparse.ArgumentParser(description='Program to take \
        dictionaries of variants/consistency counts and a coverage threshold \
        and returns a dataframe with all types of consistency (to add: \
        filter by other criteria, such as imputation R^2/coverage).')

parser.add_argument('-i', '--indir', type=str, required=True, \
        help='Input directory')

parser.add_argument('-p', '--pkldir', type=str, required=True, \
        help='Directory with consistency dicts')

parser.add_argument('-d', '--donors', type=str, required=True, \
        help='txt file containing list of multiplexed donors')

parser.add_argument('-c', '--cov_thresh', type=str, required=True, \
        help='experiment-wide coverage threshold')

parser.add_argument('-o', '--outdir', type=str, required=True, \
        help='directory to store sparse matrices per individual')

args = parser.parse_args()

indir = args.indir
pkldir = args.pkldir
donor_path = args.donors
cov_thresh = int(args.cov_thresh)
outdir = args.outdir

###### filter SNPs ######
### to do: incorporate other filters, such as R^2

def filter_snps(snps):
  pass

###### main I/O ######

donors = list(np.loadtxt(donor_path, dtype=str))
n_donors = len(donors)
barcodes = pd.Index(pd.read_csv(f'{indir}/barcodes.tsv.gz',
                       sep='\t',header=None,index_col=0).index)

dp = mmread(f'{indir}/cellSNP.tag.DP.mtx').tocsr()

vcf = pd.read_csv(f'{indir}/varcon.SNPs.vcf.gz', sep='\t', header=0, index_col=0)
vcf['chrom_pos'] = vcf.index
vcf.reset_index(drop=True, inplace=True)
vcf['DP'] = [int(x.split(';')[1].split('=')[1]) for x in vcf['INFO']]

consistent = [mmread(f'{indir}/{donor}.consistent.mtx').tocsr()
              for donor in donors]
inconsistent = [dp - mtx for mtx in consistent]

with open(f'{pkldir}/c1_dict.pkl', 'rb') as f:
    c1_dict = pickle.load(f)

with open(f'{pkldir}/c2_dict.pkl', 'rb') as f:
    c2_dict = pickle.load(f)

with open(f'{pkldir}/i1_dict.pkl', 'rb') as f:
    i1_dict = pickle.load(f)

with open(f'{pkldir}/i2_dict.pkl', 'rb') as f:
    i2_dict = pickle.load(f)

###

bcs_donors = [[f'{bc}_{donor}' for donor in donors] for bc in barcodes]
c1_indices = [[c1_dict[x] for x in y] for y in bcs_donors]
c2_indices = [[c2_dict[x] for x in y] for y in bcs_donors]
i1_indices = [[i1_dict[x] for x in y] for y in bcs_donors]
i2_indices = [[i2_dict[x] for x in y] for y in bcs_donors]

### opportunity to filter variants

n_bcs = barcodes.shape[0]

c1_counts = []
c2_counts = []
i1_counts = []
i2_counts = []

for i in tqdm(range(n_bcs), total=n_bcs):
    tmp_c1 = c1_indices[i]
    tmp_c2 = c2_indices[i]
    tmp_i1 = i1_indices[i]
    tmp_i2 = i2_indices[i]

    tmp_c1_counts = []
    tmp_c2_counts = []
    tmp_i1_counts = []
    tmp_i2_counts = []

    for j in range(n_donors):
        ### do some sort of variant filtering here
        c1_mask = vcf.iloc[tmp_c1[j]]['DP'] > cov_thresh
        c2_mask = vcf.iloc[tmp_c2[j]]['DP'] > cov_thresh
        i1_mask = vcf.iloc[tmp_i1[j]]['DP'] > cov_thresh
        i2_mask = vcf.iloc[tmp_i2[j]]['DP'] > cov_thresh

        final_c1_mask = vcf.iloc[tmp_c1[j]][c1_mask].index
        final_c2_mask = vcf.iloc[tmp_c2[j]][c2_mask].index
        final_i1_mask = vcf.iloc[tmp_i1[j]][i1_mask].index
        final_i2_mask = vcf.iloc[tmp_i2[j]][i2_mask].index

        tmp_c1_counts.append(np.sum(consistent[j][final_c1_mask, i]))
        tmp_c2_counts.append(np.sum(consistent[j][final_c2_mask, i]))
        tmp_i1_counts.append(np.sum(inconsistent[j][final_i1_mask, i]))
        tmp_i2_counts.append(np.sum(inconsistent[j][final_i2_mask, i]))

    c1_counts.append(tmp_c1_counts)
    c2_counts.append(tmp_c2_counts)
    i1_counts.append(tmp_i1_counts)
    i2_counts.append(tmp_i2_counts)

c1_df = pd.DataFrame(c1_counts, columns=donors, index=barcodes)
c2_df = pd.DataFrame(c2_counts, columns=donors, index=barcodes)
i1_df = pd.DataFrame(i1_counts, columns=donors, index=barcodes)
i2_df = pd.DataFrame(i2_counts, columns=donors, index=barcodes)

###

cov = str(cov_thresh)
c1_df.to_csv(f'{outdir}/c1_df_{cov}.csv', sep='\t', header=True, index=True)
c2_df.to_csv(f'{outdir}/c2_df_{cov}.csv', sep='\t', header=True, index=True)
i1_df.to_csv(f'{outdir}/i1_df_{cov}.csv', sep='\t', header=True, index=True)
i2_df.to_csv(f'{outdir}/i2_df_{cov}.csv', sep='\t', header=True, index=True)
