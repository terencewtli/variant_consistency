#!/usr/bin/env python
import sys
import pandas as pd

### Summarise the toy example: per-droplet I1 rate for the toy method's singlet calls, split by whether the call was correct.
### usage: python summarize_toy.py <per-method table, e.g. toy/varcon/per_method/toy_method_0.csv> <toy/toy_method_status.tsv>

varcon = pd.read_csv(sys.argv[1], sep='\t', index_col=0)
status = pd.read_csv(sys.argv[2], sep='\t', index_col=0)['status']
cats = ['C1', 'C2', 'I1', 'I2']
varcon = varcon[varcon[cats].sum(axis=1) > 0]
rates = varcon[cats].div(varcon[cats].sum(axis=1), axis=0)
rates['status'] = status.reindex(rates.index)
summary = rates.groupby('status')[cats].median().mul(100).round(1)
summary.insert(0, 'droplets', rates.groupby('status').size())
print('median % of allele counts per droplet, toy_method singlet calls:')
print(summary.to_string())
