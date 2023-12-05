# THIS PROGRAM INCLUDES TWO MAIN FUNCTIONS:
#    1) WRITE_READS_OUTPUT(): TO DETERMINE THE READS GIVEN A SET OF UORF COORDINATES
#    2) WRITE_PERIODICITY_OUTPUT(): TO DETERMINE THE PERIODICITY OF A SET OF UORFS

import ribopy
import ribopy.rnaseq
import numpy as np
import pandas as pd
from ribopy import Ribo
import scipy
import sys
import csv
import copy
import matplotlib.pyplot as plt
import gzip
from scipy.stats import chisquare
from ast import literal_eval
from statsmodels.stats.multitest import multipletests

#### to add a key-value pair to a dictionary
def add_dic(key, value, dict):
    if key in dict.keys():
        dict[key].append(value)
    else:
        dict[key] = [value]

#### function to determine periodicity output using chi-squared and FDR
def find_periodicity(df):
    frames_columns = df.filter(like="Frames", axis=1)
    frames_column_names = frames_columns.columns.tolist()

    # Create a long-format DataFrame using melt
    long_df = pd.melt(frames_columns.reset_index(), id_vars='index', value_vars=frames_column_names, var_name="Experiment", value_name="Frames")

    # Remove "Frames " from the Experiment column
    long_df["Experiment"] = long_df["Experiment"].str.replace("Frames", "")

    # Convert the string representations of lists to actual lists
    long_df["Frames"] = long_df["Frames"].apply(literal_eval)

    # Calculate chi-square p-values with FDR
    output_sample = {**long_df, **find_chisquared(long_df["Frames"])}

    # Data organization
    sample_df = pd.DataFrame({
        "Index": output_sample["index"],
        "Experiment": output_sample["Experiment"],
        "Reject": output_sample["Reject"],
        "Corrected p-value": output_sample["Corrected p-value"]
    })

    df_reject = sample_df.pivot(index=["Index"], columns='Experiment', values='Reject')
    df_p_value = sample_df.pivot(index=["Index"], columns='Experiment', values='Corrected p-value')

    # Rename the columns
    df_reject.columns = [f"{col} Reject" for col in df_reject.columns]
    df_p_value.columns = [f"{col} corrected p-value" for col in df_p_value.columns]

    # Concatenate the two DataFrames
    final_df = pd.concat([df_reject, df_p_value], axis=1)
    final_df = pd.concat([df, final_df], ignore_index=False, axis=1)

    return final_df

#### performs chi-squared with FDR correction
def find_chisquared(arr):
    p_vals = []
    output = {}
    for row in arr:
        total = sum(row)
        chi2, p_value = chisquare(row, f_exp=[total / 3, total / 3, total / 3])
        p_vals.append(p_value)

    p_vals = np.array(p_vals)
    mask = np.isfinite(p_vals)
    pval_corrected = np.empty(p_vals.shape)
    pval_corrected.fill(np.nan)

    reject_list = np.empty(p_vals.shape)
    reject_list.fill(np.nan)

    reject_list[mask], pval_corrected[mask] = multipletests(p_vals[mask], method='fdr_bh', alpha=0.05)[:2]

    output["Reject"] = reject_list.tolist()
    output["Corrected p-value"] = pval_corrected

    return output

#### calls relevant functions to output reads
# input:
#    transcripts = csv of uORF to generate reads for, e.g. "ltdstart_uorf_no_overlap.csv"
#    outfile = name for outputted csv file
def write_reads_output(transcripts, out_file):
    transcript_df = pd.read_csv(transcripts)


    reads = find_reads(transcript_df)
    reads = pd.DataFrame(reads)

    output = pd.concat([transcript_df, reads], ignore_index=False, axis=1)
    output.to_csv(out_file, index=False)

# calls relevant functions to output periodicity results
# inputs:
#   file = csv of filtered uORFs, must have a Frames column
#   complete_outfile = name of output csv file with all periodicity results
#   accepted_outfile = name of output csv file with periodicity results for uORFs determined to be periodic
def write_periodicity_output(file, complete_outfile, accepted_outfile):
    df = pd.read_csv(file, sep=',')
    df2 = find_periodicity(df)
    df2.to_csv(complete_outfile, index=False)

    # Filter rows where the sum of values in "Reject" columns is greater than 0
    columns_containing_reject = [col for col in df2.columns if 'Reject' in col]

    filtered_df = df2[df2[columns_containing_reject].sum(axis=1) > 0]
    filtered_df.to_csv(accepted_outfile, index = False)

#### TO RUN:

# to check PERIODICITY:
write_periodicity_output("UORF_LIST.csv", "COMPLETE_OUTPUT.csv", "ACCEPTED_OUTPUT.csv")


