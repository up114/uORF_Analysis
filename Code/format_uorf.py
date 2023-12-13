# Uma Paul - 1/19/23
# This program counts to number of reads for an array of genes

import ribopy
import ribopy.rnaseq
import numpy as np
import pandas as pd
from ribopy import Ribo
import sys
import csv
import copy
import gzip
from scipy.stats import chisquare
from statsmodels.stats.multitest import multipletests

ribo_path = "merged_mapq10.ribo"

### code to find CDS coordinates
def find_CDS(genes) :
    output = {}
    transcr_names = {}
    ribo_object = Ribo(ribo_path)
    transcripts = ribo_object.transcript_names
    for transcript in transcripts :
        split_transcript = transcript.split("|")
        gene = split_transcript[4]
        if gene in genes:
            CDS = split_transcript[8][4:].split("-")
            output[gene] = [int(CDS[0]), int(CDS[1])]
            trans_name = split_transcript[0]
            transcr_names[gene] = trans_name
    return (output, transcr_names)

# add value-key pairs to a dictionary
def add_dic(key, value, dict):
    if key in dict.keys():
        dict[key].append(value)
    else:
        dict[key] = [value]

# codons to amino acids
def translate_sequence(sequence):
    genetic_code = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    amino_acids = []
    codons = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]

    for codon in codons:
        amino_acid = genetic_code.get(codon.upper(), '')
        amino_acids.append(amino_acid)

    return ''.join(amino_acids)

# adds relevant information
def write_lab_output(transcripts, outfile) :
    output = {}

    with gzip.open('appris_mouse_v2_selected.fa.gz', 'rb') as f:
        sequence = str(f.read())

    transcript_df = pd.read_csv(transcripts)
    genes = transcript_df["Gene"].values.tolist()
    starts = transcript_df["Start"].values.tolist()
    stops = transcript_df["Actual Stop"].values.tolist()

    cds, trans_names = find_CDS(genes)

    for i in range(len(genes)) :
        gene = genes[i]
        start = starts[i]
        stop = stops[i]
        loc = sequence.find(gene)
        file_name = sequence[loc:sequence.find("\\n", loc)]

        current = sequence[loc + len(file_name): sequence.find(">", loc + 1)].replace("\\n", "")
        nt = current[start : stop]
        comp_seq = current[0:stop - 3]
        add_dic("Transcript", trans_names[gene], output)
        add_dic("Kozak", current[start - 3:start + 4], output)
        add_dic("Sequence", nt, output)
        add_dic("Complete Sequence", comp_seq, output)
        add_dic("AA Sequence", translate_sequence(nt), output)
        add_dic("Length", stop-start, output)
        add_dic("AA Length", (stop-start) / 3, output)
        add_dic("CDS Start", cds[gene][0], output)
        add_dic("CDS Stop", cds[gene][1], output)

    out_file = pd.concat([pd.read_csv(transcripts), pd.DataFrame.from_dict(output)], axis = 1)
    out_file.to_csv(outfile, index = False)

# TO RUN:
write_lab_output("INPUT.csv", "OUTPUT.csv")
