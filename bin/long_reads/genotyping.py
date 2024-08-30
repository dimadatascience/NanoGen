#!/usr/bin/python

import sys
import pandas as pd
import numpy as np
from scipy.stats import nbinom


##


# Main
def main():

    path_input = sys.argv[1]
    sample_name = sys.argv[2]
    # path_input = '/Users/IEO5505/Desktop/scmseq_new/scratch/counts.tsv.gz'

    pseudocount = .000000001
    barcodes = pd.read_csv(path_input).rename(columns={'barcode':'cell'})
    grouped = barcodes.groupby('cell')
    barcodes = barcodes.merge(    
        ( 
            (grouped['MIS'].sum() / 
            (grouped[['MUT', "WT", "MIS"]].sum().sum(axis=1) + pseudocount) / 2)
            .to_frame('err').reset_index()
        ),
        on='cell'
    )
    barcodes.loc[barcodes['err']==0, 'err'] = barcodes['err'].min()

    mut = barcodes["MUT"].values
    wt = barcodes["WT"].values
    experimental_error = barcodes["err"].values
    experimental_error[experimental_error<1e-4] = 1e-4
    p = nbinom.cdf(k=wt, n=mut, p=experimental_error)
    barcodes['p'] = p

    genotypes = np.zeros(p.size)
    genotypes[(wt==0) & (mut==0)] = np.nan      # Real NAs: counts mut==0 and counts wt==0
    genotypes[p<.1] = 1                         # MUT: p WT very low (<.1)
    genotypes[(np.isnan(p)) & (wt>=3)] = 0      # WT: p NA for counts mut == 0 and count wt>=3
    genotypes[p>=.9] = 0                        # WT: p WT very high (>=.9)
    genotypes[(np.isnan(p)) & (wt<3)] = np.nan  # NAs (impossible to assign): p NA for counts mut == 0 but count wt<3
    genotypes[(p>=.1) & (p<.9)] = np.nan        # NAs (impossible to assign): p WT between .1 and .9

    barcodes['genotype'] = np.select([genotypes==0, genotypes==1], ['WT', 'MUT'], default='-')
    barcodes.loc[barcodes['genotype']=='-','genotype'] = np.nan
    barcodes.to_csv(f'{sample_name}.csv')


##


# Run
if __name__ == "__main__":
    main()
