# THIS CODE IS USED TO COMPUTATIONALLY GENERATE ALL POSSIBLE UORF CANDIDATES BASED ON THE APPRIS TRANSCRIPT
import gzip # to unzip the gz file
import re
import numpy as np
import pandas as pd

### code to find CDS start and stop
def find_CDS(transcript_name) :
    split_transcript = transcript_name.split("|")
    CDS = split_transcript[8][4:].split("-")
    return split_transcript[4], int(CDS[0]), int(CDS[1])


# add key-value pair to dictionary
def add_dic(key, value, dict):
    if key in dict.keys():
        dict[key].append(value)
    else:
        dict[key] = [value]

# find indexes of substrings in a string
def find_multiple_substrings(string, substrings):
    pattern = '|'.join(map(re.escape, substrings))
    pattern = f'(?=({pattern}))'  # Add a positive lookahead assertion
    matches = re.finditer(pattern, string)
    return [match.start(1) for match in matches]

# find all uORFs in APPRIS
def find_uORFS( min, max, sequence ) :
    uorf_list = {} # includes complete uORFs with CDS overlapping portions
    uorf_list_no_overlap = {} # uORFs with CDS overlap removed

    file_start = sequence.find(">",0)
    curr = 1
    while (file_start != -1):
        curr += 1
        next_file_start = sequence.find(">", file_start + 1)
        file_name = sequence[file_start + 1 :sequence.find("\\n", file_start)]

        if(file_name.find('UTR5') == -1) :
            file_start = next_file_start
            continue

        gene, CDS_start, CDS_stop = find_CDS(file_name)

        current = sequence[file_start + len(file_name) + 1:next_file_start].replace("\\n", "")
        starts = np.array(find_multiple_substrings(current, ["ATG", "ACG", "ATT", "CTG", "GTG", "TTG"]))
        stops = np.array(find_multiple_substrings(current, ["TAA", "TGA", "TAG"]))
        starts = starts[starts < (CDS_start - 1)]
        stops = stops[stops + 3 <= CDS_stop]
        for start in starts:
            for stop in stops[stops > start] :
                stop = stop + 3
                orf_len = stop - start
                if(orf_len % 3 == 0) :
                    if (orf_len >= min and orf_len <= max) :
                        add_dic("Transcript", file_name.split("|")[0], uorf_list)
                        add_dic("Gene", gene, uorf_list)
                        add_dic("Start", start, uorf_list)
                        add_dic("Stop", stop, uorf_list)
                        add_dic("Start Codon", current[start:start + 3], uorf_list)
                        add_dic("CDS stop", CDS_stop, uorf_list)

                        if stop > CDS_start :
                            if stop == CDS_stop:
                                add_dic("Gene", gene, uorf_list_no_overlap)
                                add_dic("Start", start, uorf_list_no_overlap)
                                add_dic("Stop", CDS_start - 1, uorf_list_no_overlap)
                                add_dic("Start Codon", current[start:start + 3], uorf_list_no_overlap)
                                add_dic("Actual Stop", stop, uorf_list_no_overlap)
                                add_dic("Type", "extension", uorf_list_no_overlap)

                            else:
                                add_dic("Gene", gene, uorf_list_no_overlap)
                                add_dic("Start", start, uorf_list_no_overlap)
                                add_dic("Stop", CDS_start - 1, uorf_list_no_overlap)
                                add_dic("Start Codon", current[start:start + 3], uorf_list_no_overlap)
                                add_dic("Actual Stop", stop, uorf_list_no_overlap)
                                add_dic("Type", "overlapping", uorf_list_no_overlap)
                        else :
                            add_dic("Gene", gene, uorf_list_no_overlap)
                            add_dic("Start", start, uorf_list_no_overlap)
                            add_dic("Stop", stop, uorf_list_no_overlap)
                            add_dic("Start Codon", current[start:start + 3], uorf_list_no_overlap)
                            add_dic("Actual Stop", stop, uorf_list_no_overlap)
                            add_dic("Type", "nonoverlapping", uorf_list_no_overlap)
                    break
        file_start = next_file_start
    return uorf_list, uorf_list_no_overlap

# find all uORFs and generate a csv file with uORF sequence with and without overlap
def write_uORFs(file, min_length, max_length) :
    # read the file
    with gzip.open('appris_mouse_v2_selected.fa.gz', 'rb') as f:
        sequence = str(f.read())  # contains a string with the entire listing of files

    uorf, uorf_nooverlap = find_uORFS(min_length, max_length, sequence)
    uorf_df = pd.DataFrame.from_dict(uorf)
    uorf_overlap_df = pd.DataFrame.from_dict(uorf_nooverlap)
 #   uorf_df.to_csv("ltdstart_uorf_transcript_all.csv", index=False)
    uorf_overlap_df.to_csv(file, index = False)

# TO RUN:
write_uORFs("FILE.csv", 33, 300)








