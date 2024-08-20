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
    geno_pred = f'{sample_name}.xlsx'
    tsv_pred = f'{sample_name}.tsv'

    barcodes = pd.read_csv(path_input).rename(columns={'barcode':'cell'})
    barcodes = (
        barcodes.merge(
            barcodes.groupby('cell')
            .apply(lambda x: x['MIS'].values.sum() / x[['MUT', "WT", "MIS"]].values.sum() / 2).to_frame('err').reset_index(),
            on='cell'
        )
    )
    barcodes.loc[barcodes['err']==0, 'err'] = barcodes['err'].min()

    mut = barcodes["MUT"].values
    wt = barcodes["WT"].values
    experimental_error = barcodes["err"].values
    p = nbinom.cdf(k=wt, n=mut, p=experimental_error)
    barcodes['p'] = p

    genotypes = np.zeros(p.size)
    genotypes[(wt==0) & (mut==0)] = np.nan      # Real NAs: counts mut==0 and counts wt==0
    genotypes[p<.1] = 1                         # MUT: p WT very low (<.1)
    genotypes[(np.isnan(p)) & (wt>=3)] = 0      # WT: p NA for counts mut == 0 and count wt>=3
    genotypes[p>=.9] = 0                        # WT: p WT very high (>=.9)
    genotypes[(np.isnan(p)) & (wt<3)] = np.nan  # NAs (impossible to assign): p NA for counts mut == 0 but count wt<3
    genotypes[(p>=.1) & (p<.9)] = np.nan        # NAs (impossible to assign): p WT between .1 and .9

    barcodes['genotype'] = np.select([genotypes==0, genotypes==1], ['WT', 'MUT'], default=np.nan)
    barcodes.to_excel(geno_pred)
    barcodes.to_csv(tsv_pred)


##


# Run
if __name__ == "__main__":
    main()
