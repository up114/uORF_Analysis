# THIS PROGRAM COUNTS THE READS GIVEN A SET OF UORFS:
#    1) WRITE_READS_OUTPUT(): TO DETERMINE THE READS GIVEN A SET OF UORF COORDINATES

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
        max_value_index = row.iloc[35:41].idxmax() # p-site offsets are restricted to an 11-16 nt range
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

#### TO RUN:
# set ribo file
ribo_path = "neural.ribo"

# to check READS:
write_reads_output("UORFS_LIST.csv", "OUTPUT.csv")



