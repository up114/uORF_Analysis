# This program prints the longest open reading frame (ORF) for a file containing multiple mRNA sequences,
# each with a unique file name

import gzip # to unzip the gz file
import re
import numpy as np
import pandas as pd

### code to find CDS start and stop
def find_CDS(transcript_name) :
    split_transcript = transcript_name.split("|")
    CDS = split_transcript[8][4:].split("-")
    return split_transcript[4], int(CDS[0]), int(CDS[1])

### remove CDS overlap from starts/stops
def remove_CDS(genes, starts, stops) :
    CDS = find_CDS(genes)
    out_start = starts[:]
    out_stop = stops[:]
    print(genes)
    print(starts)
    print(stops)
    for i in range(len(genes)):
        gene_CDS = CDS[genes[i]][0]
        if(starts[i] > gene_CDS[0] - 15 and starts[i] < gene_CDS[1] + 15):
            out_start[i] = gene_CDS[1] + 15
        if (stops[i] > gene_CDS[0] - 15 and stops[i] < gene_CDS[1] + 15):
            out_stop[i] = gene_CDS[0] - 15
    return out_start, out_stop

def add_dic(key, value, dict):
    if key in dict.keys():
        dict[key].append(value)
    else:
        dict[key] = [value]

def find_multiple_substrings(string, substrings):
    pattern = '|'.join(map(re.escape, substrings))
    matches = re.finditer(pattern, string)
    return [match.start() for match in matches]

def find_uORFS( min, max, sequence ) :
    uorf_list = {} # includes complete uORFs with CDS overlapping portions
    uorf_list_no_overlap = {} # uORFs with CDS overlap removed

    file_start = sequence.find(">",0)

    while (file_start != -1):
        next_file_start = sequence.find(">", file_start + 1)
        file_name = sequence[file_start:sequence.find("\\n", file_start)]

        if(file_name.find('UTR5') == -1) :
            file_start = next_file_start
            continue

        gene, CDS_start, CDS_stop = find_CDS(file_name)
        print(file_name)
        current = sequence[file_start + len(file_name):next_file_start].replace("\\n", "")
        starts = np.array(find_multiple_substrings(current, ["ATG", "AAG", "ACG", "AGG", "ATA", "ATC", "ATT", "CTG", "GTG", "TTG"]))
        stops = np.array(find_multiple_substrings(current, ["TAA", "TGA", "TAG"]))
        starts = starts[starts < CDS_start]
        stops = stops[stops < CDS_stop]

        for start in starts:
            stop_passed = False
            for stop in stops[stops > start] :
                orf_len = stop - start
                if(orf_len % 3 == 0 and orf_len >= min and orf_len <= max) :
                    add_dic("Gene", gene, uorf_list)
                    add_dic("Start", start, uorf_list)
                    add_dic("Stop", stop, uorf_list)
                    if stop > CDS_start :
                        if(stop_passed == False) :
                            add_dic("Gene", gene, uorf_list_no_overlap)
                            add_dic("Start", start, uorf_list_no_overlap)
                            add_dic("Stop", CDS_start, uorf_list_no_overlap)
                            stop_passed = True
                    else :
                        add_dic("Gene", gene, uorf_list_no_overlap)
                        add_dic("Start", start, uorf_list_no_overlap)
                        add_dic("Stop", stop, uorf_list_no_overlap)

        file_start = next_file_start
    return uorf_list, uorf_list_no_overlap

#     for start in starts:
#         j = start  # stores the location of the current codon
#         while (j + 3 <= CDS_start) :
#             codon = current[j:j+3]
#             if (codon == "TGA" or codon == "TAA" or codon == "TAG"): # continue parsing until a stop codon is reached
#                 ### run 3 nt periodicity to determine eligibility ###
#                 orfs.append((start, j + 2)) ### add to uORFs if the ORF shows 3 nt periodicity
#
#             j += 3
#







# # for each mRNA sequence, find all uORFs
# for file in range(len(start_indexes) - 1):
#     # find and print the file name
#     file_name = sequence[start_indexes[file]:sequence.find("\\n", start_indexes[file])]
#     print(file_name)
#
#     # save the current mRNA sequence in current
#     current = sequence[start_indexes[file] + len(file_name):start_indexes[file + 1]]
#     current = current.replace("\\n", "")
#
#     # find all of the start codons in current and store them in start_files
#     index = 0 # stores current index while parsing
#     start_codons = [] # stores index of all start codons
#     while (index < len(current)):
#         found = current.find("ATG", index) # find the next start codon after the current index
#         if (found != -1): # if a start codon is found, save its location and increment index
#             start_codons.append(int(found))
#             index = found + 1
#         else: # if no more start codons can be found, exit the loop
#             break
#
#     longest = [0,0] # stores the longest ORF
#     orfs = []
#     end_reached = False # stores whether the end of the mRNA sequence has been reached
#     for i in range(0, len(start_codons) - 1): # parse through current for each start codon
#         codon = current[start_codons[i]:start_codons[i] + 3] # stores the current codon
#         j = start_codons[i] # stores the location of the current codon
#         while (not (codon == "TGA" or codon == "TAA" or codon == "TAG")): # continue parsing until a stop codon is reached
#             if(j >= len(current)) : # if the end of the file is reached, exit the loop
#                 end_reached = True
#                 break
#
#             # continue moving through the current sequence by incrementing j and updating codon
#             j += 3
#             codon = current[j: j + 3]
#
#         if(not (end_reached)) : # only update longest if the end of the file was not reached
#             # check if the current ORF is longer than the longest ORF stored in longest
#             orf = current[start_codons[i]: (j + 3)]
#             orfs.append((start_codons[i]+1, j+3))
#             if (len(orf) > (longest[1] - longest[0])):
#                 longest = (start_codons[i] + 1, (j + 3))
#
#     # print the longest ORF for each file
#     print("The longest orf: ", longest)
#     if(longest != [0,0]) :
#         print("The preceding orfs: ", orfs[:(orfs.index(longest))])

def write_uORFs() :
    # read the file
    with gzip.open('appris_mouse_v2_selected.fa.gz', 'rb') as f:
        sequence = str(f.read())  # contains a string with the entire listing of files

    uorf, uorf_nooverlap = find_uORFS(6, 300, sequence)
    uorf_df = pd.DataFrame.from_dict(uorf)
    uorf_overlap_df = pd.DataFrame.from_dict(uorf_nooverlap)
    uorf_df.to_csv("uorfs_lab.csv", index=False)
    uorf_overlap_df.to_csv("uorf_no_overlap_lab.csv", index = False)

write_uORFs()







