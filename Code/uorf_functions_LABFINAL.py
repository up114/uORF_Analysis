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

####create Ribo object from path
def ribo_data(ribo_path):
    ribo_object = Ribo(ribo_path, alias=defalternative_human_alias)
    return ribo_object

#### find RPF length range
def intevl(experiment_path, experiment_id):
    ribo_object = ribo_data(experiment_path)
    data_tmp = ribo_object.get_length_dist("CDS")
    data_tmp.reset_index(inplace=True)
    data = data_tmp.iloc[6:26]
    pct_85 = sum(data["%s" % experiment_id]) * 0.85
    # pct_90=sum(data["%s"%experiment_id])*0.90
    value = data[data["%s" % experiment_id] == data["%s" % experiment_id].max()]["%s" % experiment_id].values[0]
    mmin = mmax = data[data["%s" % experiment_id] == data["%s" % experiment_id].max()]['read_length'].values[0]
    while value <= pct_85:
        if mmax < 40 and mmin > 21:
            if data[data['read_length'] == mmax + 1]["%s" % experiment_id].values[0] >= \
                    data[data['read_length'] == mmin - 1]["%s" % experiment_id].values[0]:
                mmax += 1
                value += data[data['read_length'] == mmax]["%s" % experiment_id].values[0]
            else:
                mmin -= 1
                value += data[data['read_length'] == mmin]["%s" % experiment_id].values[0]
        elif mmax == 40:
            mmin -= 1
            value += data[data['read_length'] == mmin]["%s" % experiment_id].values[0]
        elif mmin == 21:
            mmax += 1
            value += data[data['read_length'] == mmax]["%s" % experiment_id].values[0]
    # print(min,max)
    read_pct = value / sum(data["%s" % experiment_id])  # to check if study is suitable for analysis
    return mmin, mmax, read_pct


#### function to create aliases and shorten transcript names
def defalternative_human_alias(x):
    x_pieces = x.split("|")
    return x_pieces[4]

#### to add a key-value pair to a dictionary
def add_dic(key, value, dict):
    if key in dict.keys():
        dict[key].append(value)
    else:
        dict[key] = [value]


####return p-site offset for an experiment
def psite_offset(ribo_object, exp, mmin, mmax):
    df = (ribo_object.get_metagene("start", experiments=exp, range_lower=mmin, range_upper=mmax, sum_lengths=False,
                                   sum_references=True))

    p_site = {}

    for index, row in df.iterrows():
        max_value_index = row.iloc[35:41].idxmax() # p-site offsets are restricted to a 11-16 nt range
        print(row.iloc[35:41])
        offset = -1 * max_value_index + 1
        p_site[index[1]] = offset

    return p_site

#### function to find reads given genes, start + stop coordinates
# input:
#   transcript_df: pandas df with "Gene", "Start", and "Stop" columns
def find_reads(transcript_df):
    genes = transcript_df["Gene"].values.tolist()
    starts = transcript_df["Start"].values.tolist()
    stops = transcript_df["Stop"].values.tolist()

    output = {}
    ribo = Ribo(ribo_path, alias=defalternative_human_alias)

    ribo.print_info()

    for experiment in ribo.experiments:

        mmin, mmax, read_pct = intevl(ribo_path, experiment)
        total = ribo.info['experiments'][experiment]['Reads']
        offset = psite_offset(ribo, experiment, mmin, mmax)

        print(f"{experiment} : {offset}")
        uorf = {}
        frame1_reads = {}
        frame2_reads = {}
        frame3_reads = {}

        for k in range(mmin, mmax + 1):
            df = ribo.get_coverage(experiment=experiment, range_lower=k, range_upper=k, alias=True)

            for l in range(len(genes)):
                gene = genes[l]
                start = starts[l]
                stop = stops[l]

                if start <= offset[k]:
                    if k == mmax:
                        add_dic('Experiment' + experiment, None, output)
                        add_dic('Frames' + experiment, None, output)
                    continue

                try:
                    coverage = df[gene]
                    reads = np.sum(coverage[start - offset[k]: stop - offset[k]])

                    if (gene, start, stop) in uorf.keys():
                        uorf[(gene, start, stop)] += reads
                        frame1_reads[(gene, start, stop)] += np.sum(coverage[start - offset[k]: stop - offset[k]: 3])
                        frame2_reads[(gene, start, stop)] += np.sum(
                            coverage[start - offset[k] + 1: stop - offset[k]: 3])
                        frame3_reads[(gene, start, stop)] += np.sum(
                            coverage[start - offset[k] + 2: stop - offset[k]: 3])
                    else:
                        uorf[(gene, start, stop)] = reads
                        frame1_reads[(gene, start, stop)] = np.sum(
                            coverage[start - offset[k]: stop - offset[k]: 3])
                        frame2_reads[(gene, start, stop)] = np.sum(
                            coverage[start - offset[k] + 1: stop - offset[k]: 3])
                        frame3_reads[(gene, start, stop)] = np.sum(
                            coverage[start - offset[k] + 2: stop - offset[k]: 3])

                    if k == mmax:
                        add_dic('Experiment' + experiment, uorf[(gene, start, stop)] * 1000000 / total, output)
                        add_dic('Frames' + experiment, [frame1_reads[(gene, start, stop)],
                                         frame2_reads[(gene, start, stop)],
                                         frame3_reads[(gene, start, stop)]], output)
                except KeyError:
                    continue
    return output

#### function to determine periodicty output using chi-sqaured and FDR
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
#    transcripts = csv of uORF to generate reads for, "ltdstart_uorf_no_overlap.csv"
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
#   outfile = name of outputted csv file
def write_periodicity_output(file, complete_outfile, accepted_outfile):
    df = pd.read_csv(file, sep=',')
    df2 = find_periodicity(df)
    df2.to_csv(complete_outfile, index=False)

    # Filter rows where the sum of values in "Reject" columns is greater than 0
    columns_containing_reject = [col for col in df2.columns if 'Reject' in col]

    filtered_df = df2[df2[columns_containing_reject].sum(axis=1) > 0]
    filtered_df.to_csv(accepted_outfile, index = False)

#### TO RUN:
# set ribo file
ribo_path = "neural.ribo"

# to check READS:
write_reads_output("ltdstart_uorf_no_overlap.csv", "neural_output_NEW.csv")

# to check PERIODICITY:
#write_periodicity_output("FILTER_OUTPUT.csv", "PERIODICITY_OUTPUT.csv")
#write_periodicity_output("neural_filterNEW.csv", "neural_periodicityTEST.csv", "neural_acceptedTEST.csv")


