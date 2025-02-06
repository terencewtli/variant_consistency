#!/usr/bin/env python

import numpy as np
import pandas as pd

import argparse
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

parser = argparse.ArgumentParser(description='Program \
        to take as input consistency counts and a demultiplexing \
        method dataframe, and returns a per-method consistency counts \
        dataframe.')

parser.add_argument('-i', '--indir', type=str, required=True, \
        help='Input directory')

parser.add_argument('-d', '--donors', type=str, required=True, \
        help='txt file containing list of multiplexed donors')

parser.add_argument('-c', '--cov_thresh', type=str, required=True, \
        help='experiment-wide coverage threshold')

parser.add_argument('-v', '--demux', type=str, required=True, \
        help='csv file containing barcodes and their demux \
        method assignment')

parser.add_argument('-o', '--outdir', type=str, required=True, \
        help='directory to store final per-method consistency')

args = parser.parse_args()

indir = args.indir
donor_path = args.donors
cov_thresh = int(args.cov_thresh)
demux_path = args.demux
outdir = args.outdir

###### main I/O ######

donors = list(np.loadtxt(donor_path, dtype=str))
n_donors = len(donors)

c1_df = pd.read_csv(f'{indir}/c1_df_{cov_thresh}.csv', sep='\t', header=0, index_col=0)
c2_df = pd.read_csv(f'{indir}/c2_df_{cov_thresh}.csv', sep='\t', header=0, index_col=0)
i1_df = pd.read_csv(f'{indir}/i1_df_{cov_thresh}.csv', sep='\t', header=0, index_col=0)
i2_df = pd.read_csv(f'{indir}/i2_df_{cov_thresh}.csv', sep='\t', header=0, index_col=0)

demux_df = pd.read_csv(demux_path, sep='\t', header=0, index_col=0)

singlet_df = demux_df.copy()
for donor in donors:
        singlet_df.replace({donor : 1}, inplace=True)
singlet_df.replace({'doublet' : 0, 'unassigned' : 0}, inplace=True)

##### get per-method counts for all barcodes

cols = ['C1', 'C2', 'I1', 'I2', 'donor']
demux_methods = demux_df.columns

for method in demux_methods:
    nonsing_mask = singlet_df[method] == 1
    singlet_bcs = singlet_df[nonsing_mask].index
    varcon_df = pd.DataFrame(index = singlet_bcs, columns = cols)
    for barcode in singlet_bcs:
        donor = demux_df.loc[barcode, method]
        c1_counts = c1_df.loc[barcode][donor]
        c2_counts = c2_df.loc[barcode][donor]
        i1_counts = i1_df.loc[barcode][donor]
        i2_counts = i2_df.loc[barcode][donor]
        varcon_df.loc[barcode] = [c1_counts, c2_counts, i1_counts, i2_counts, donor]
    varcon_df.to_csv(f'{outdir}/{method}_{cov_thresh}.csv', sep='\t', header=True, index=True)

##### summarize mean rates across all barcodes and write into one df

cols = ['C1', 'C2', 'I1', 'I2']
demux_methods = demux_df.columns

demux_rates = pd.DataFrame(index=demux_methods, columns=cols)

for method in demux_methods:
    nonsing_mask = singlet_df[method] == 1
    singlet_bcs = singlet_df[nonsing_mask].index
    varcon_demux_df = pd.DataFrame(index = singlet_bcs, columns = cols)
    for barcode in singlet_bcs:
        donor = demux_df.loc[barcode, method]
        c1_counts = c1_df.loc[barcode][donor]
        c2_counts = c2_df.loc[barcode][donor]
        i1_counts = i1_df.loc[barcode][donor]
        i2_counts = i2_df.loc[barcode][donor]
        varcon_demux_df.loc[barcode] = [c1_counts, c2_counts, i1_counts, i2_counts]
    total_counts = varcon_demux_df.sum()
    denom = sum(total_counts)
    tmp_rates = total_counts.div(denom)
    demux_rates.loc[method] = tmp_rates

demux_rates.to_csv(f'{outdir}/inter_{cov_thresh}_varcon.csv', sep='\t', header=True, index=True)
