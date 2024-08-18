#!/usr/bin/python

import sys
import pandas as pd
import numpy as np
from datetime import datetime
from scipy.stats import nbinom
import matplotlib.pyplot as plt


##


# Paths and vars
mut = sys.argv[1]
sample_name = sys.argv[2]
geno_pred = f'{sample_name}.xlsx'
tsv_pred = f'{sample_name}.tsv'
log = 'genotyping.log'


##


# Logging utils
def start_log():
    with open(log, "w") as log_file:
        print("[INFO] Start the analysis", file=log_file)
        check_time()

def check_time():
    with open(log, "a") as log_file:
        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
        print("[INFO] Date and Time = {}".format(dt_string), file=log_file)


##


# Main
def main():

    # Begin logging...
    start_log()

    # TRY
    try:

        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_rows', None)

        barcodes = pd.read_csv(mut, index_col=[0])
        counts_per_cell = np.sum(barcodes.loc[:, ~(barcodes.columns.isin(['barcode']))].values, axis=1)
        counts_per_cell[counts_per_cell == 0] = np.min(counts_per_cell[counts_per_cell != 0])
        plt.hist(np.log10(counts_per_cell))
        plt.xlabel("Counts per cell log10")

        barcodes = barcodes[np.sum(barcodes.loc[:, ~(barcodes.columns.isin(['barcode']))].values, axis=1) > 3]
        genes = np.unique([x[0] for x in barcodes.columns.str.split("_").values if x[0] not in ['barcode']])

        genotypes = barcodes[['barcode']].copy()

        experimental_error = np.sum(barcodes.loc[:, ["mis" in x for x in barcodes.columns]].values, axis=1) / np.sum(
            barcodes.loc[:, ~(barcodes.columns.isin(['barcode']))].values, axis=1) / 2.
        experimental_error[experimental_error > 1] = 1
        experimental_error[experimental_error == 0] = np.min(experimental_error[experimental_error != 0])

        for idx in range(len(genes)):
            selected = barcodes.loc[:, [genes[idx] in x for x in barcodes.columns]].copy()
            genotypes[genes[idx]] = "-"

            r = selected.loc[:, ["MUT" in x for x in selected.columns]].values.squeeze()  # success
            k = selected.loc[:, ["WT" in x for x in selected.columns]].values.squeeze()  # failure

            probability_null_hypothesis_wt = nbinom.cdf(k=k, n=r, p=experimental_error)
            genotypes.loc[probability_null_hypothesis_wt < 0.1, genes[idx]] = "MUT"
            genotypes.loc[probability_null_hypothesis_wt >= 0.9, genes[idx]] = "WT"
            genotypes.loc[(np.isnan(probability_null_hypothesis_wt)), genes[
                idx]] = "WT"  # if no mutation reads is found, then nbinom is nan
            genotypes.loc[np.sum(selected, axis=1) < 3, genes[idx]] = ""

        nan_value = float("NaN")
        genotypes.replace("", nan_value, inplace=True)
        genotypes.dropna(how='all', axis=1, inplace=True)

        genotypes.to_excel(geno_pred, index=False)
        genotypes.to_csv(tsv_pred, sep="\t",index=False)

    # EXCEPT
    except Exception as e:

        with open(log, "a") as log_file:
            for line in str(e).split("\n"):
                print(f"[ERROR {line}", file=log_file)
        raise


##


# Run
if __name__ == "__main__":
    main()