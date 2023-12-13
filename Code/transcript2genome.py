#### CONVERTS TRANSCRIPT COORDINATES TO GENOMIC COORDINATES
import ribopy
import ribopy.rnaseq
import numpy as np
import pandas as pd
from ribopy import Ribo
import scipy
import sys
import csv
import copy
import gzip
import matplotlib.pyplot as plt

# add key-value pair to a dictionary
def add_dic(key, value, dict):
    if key in dict.keys():
        dict[key].append(value)
    else:
        dict[key] = [value]

# Read the GENCODE GTF file and filter for CDS entries
# input:
#    gtf_file: orginal GTF file
#    cds_output: name of output GTF file with only CDS rows
def write_cds_file(gtf_file, cds_output) :
    with open(gtf_file, 'r') as gtf_file_r:
        filtered_gtf_lines = [line for line in gtf_file_r if '\tCDS' in line]

    # Save the filtered CDS entries to a new file
    with open(cds_output, 'w') as filtered_gtf_file:
        filtered_gtf_file.writelines(filtered_gtf_lines)

# convert transcript coordinates to genomic coordinates
# inputs:
#    csv_file: file with uORF transcript coordinates; must have a "Transcript", "Start", "Stop" column
#    gtf_file: genome GTF file
#    transcript_gz: transcriptome fasta GZ file
#    outfile: name of output csv file
def transcript2genome(csv_file, gtf_file, transcript_gz, outfile) :
    # load appris file
    with gzip.open(transcript_gz, 'rb') as f:
        appris_file = str(f.read())

    output = {}

    gtf_all = pd.read_csv(gtf_file, sep='\t', header=None, comment='#')
    gtf_df = gtf_all[gtf_all[8].str.contains('appris_principal_')]
    gtf_df = gtf_df[gtf_df[2] == 'exon']

    csv_df = pd.read_csv(csv_file)
    transcripts = csv_df["Transcript"].tolist()
    starts = csv_df["Start"].tolist()
    stops = csv_df["Stop"].tolist()

    for i in range(len(transcripts)):
        transcript = transcripts[i]
        start = starts[i]
        stop = stops[i]

        gene_df = gtf_df[gtf_df[8].str.contains(transcript)]
        exon_array = []
        for _, row in gene_df.iterrows():
            exon_array.extend([i for i in range(int(row[3]), row[4] + 1)])
        exon_array.sort()
        ind = appris_file.find(transcript)
        appris_length = int(appris_file[ind: appris_file.find("\\n", ind)].split("|")[6])

        # check if array length match; if not, print the transcript and leave row empty
        if(len(exon_array) != appris_length) :
            print(transcript, len(exon_array), appris_length)
            add_dic("Chromosome", None, output)
            add_dic("Start", None, output)
            add_dic("Stop", None, output)
            add_dic("Direction", None, output)

        else:
            # for positive strands
            add_dic("Chromosome", gene_df.iloc[0, 0], output)
            if(gene_df.iloc[0, 6] == "+") :
                add_dic("Start", exon_array[int(start)], output) 
                add_dic("Stop", exon_array[int(stop)] - 1, output)
                add_dic("Direction", "+", output)
            # for negative strands
            elif(gene_df.iloc[0, 6] == "-") :
                add_dic("Start", exon_array[-int(stop) + 1] - 1, output)
                add_dic("Stop", exon_array[-int(start)] - 1, output)
                add_dic("Direction", "-", output)

    # add columns to original data and output CSV
    out_file = pd.concat([csv_df, pd.DataFrame.from_dict(output)], axis = 1)
    out_file.to_csv(outfile, index = False)

transcript2genome("ltdstart_output_lab.csv", "gencode.vM25.annotation.gtf" , "appris_mouse_v2_selected.fa.gz", "OUTPUT.csv")
