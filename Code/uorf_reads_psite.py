# THIS PROGRAM FINDS UORF READS IN RIBOBASE GIVEN A CSV WITH UORF GENES, STARTS, AND STOP (NOTE: NO HEADER IN CSV FILE)

import pandas as pd
import csv
import os
import sys
import string
import numpy as np
from ribopy import Ribo

###read csv file
def read_table(Experiment_file):
    table = pd.read_csv(Experiment_file)
    return table

####find studies and experiments to examine
def study_folder(database):
    study_dict = {}
    type_dict = {}

    study_list = read_table(database)['study_name'].values.tolist()
    experiment_list = read_table(database)['experiment_alias'].values.tolist()
    cell_line = read_table(database)['updated_cell_line_liaoyi'].values.tolist()
    for i in range(len(experiment_list)):
        if(study_list[i] != 'NaN') :
            add_dic(study_list[i], experiment_list[i], study_dict)
            type_dict[experiment_list[i]] = cell_line[i]
    return study_dict, type_dict

#####dynamic cutoff
def intevl(experiment_path, experiment_id):
    ribo_object = ribo_data(experiment_path)
    data_tmp=ribo_object.get_length_dist("CDS")
    data_tmp.reset_index(inplace=True)
    data=data_tmp.iloc[6:26]
    pct_85=sum(data["%s"%experiment_id])*0.85
    #pct_90=sum(data["%s"%experiment_id])*0.90
    value=data[data["%s"%experiment_id]==data["%s"%experiment_id].max()]["%s"%experiment_id].values[0]
    mmin=mmax=data[data["%s"%experiment_id]==data["%s"%experiment_id].max()]['read_length'].values[0]
    while value<=pct_85 :
        if mmax<40 and mmin>21:
            if data[data['read_length']==mmax+1]["%s"%experiment_id].values[0] >= data[data['read_length']==mmin-1]["%s"%experiment_id].values[0]:
                mmax+=1
                value+=data[data['read_length']==mmax]["%s"%experiment_id].values[0]
            else:
                mmin-=1
                value+=data[data['read_length']==mmin]["%s"%experiment_id].values[0]
        elif mmax==40:
            mmin-=1
            value+=data[data['read_length']==mmin]["%s"%experiment_id].values[0]
        elif mmin==21:
            mmax+=1
            value+=data[data['read_length']==mmax]["%s"%experiment_id].values[0]
    #print(min,max)
    read_pct=value/sum(data["%s"%experiment_id]) # to check if study is suitable for analysis
    return int(mmin),int(mmax),read_pct

####simplify transcript name to gene name
def defalternative_human_alias(x):
    x_pieces = x.split("|")
    return x_pieces[4]

####create Ribo object from path
def ribo_data(ribo_path):
    ribo_object = Ribo(ribo_path, alias=defalternative_human_alias)
    return ribo_object

####add a value to a dictionary
def add_dic(key, value, dict):
    if key in dict.keys():
        dict[key].append(value)
    else:
        dict[key] = [value]

####return p-site offset for an experiment
def psite_offset(ribo_object, exp, mmin, mmax) :
    df = (ribo_object.get_metagene("start", experiments = exp, range_lower= mmin, range_upper= mmax, sum_lengths = False, sum_references = True))

    p_site = {}
    #nt_per = {}
    # for index, row in df.iterrows() :
    #     max_value_index = row.iloc[35:41].idxmax()

    for index, row in df.iterrows():
        max_value_index = row.iloc[35:41].idxmax()
        offset = -1 * max_value_index + 1

        p_site[index[1]] = offset
        # ribo_object.plot_metagene("start", title = exp + " " + str(index[1]) + " : " + str(offset), experiments= exp, range_lower=index[1], range_upper=index[1],
        #                output_file= "/scratch/09369/umapaul/output/" + exp + "_" + str(index[1]) + "_start.png")
        # np_row = np.array(row)
        # frame_count = np.sum(np_row[offset % 3 :: 3])
        # total_count = np.sum(np_row)
        # nt_per[index[1]] = frame_count * 1.0 / total_count
    return p_site#, nt_per

def get_region_reads(ribo_object, region, mmin, mmax, experiment) :
    CDS_reads = (ribo_object.get_region_counts(region_name=region,
                                          range_lower=mmin,
                                          range_upper=mmax,
                                          sum_lengths=True,
                                          sum_references=False,
                                          alias=True,
                                          experiments=experiment))
    return CDS_reads


####finds the number of reads of uorfs within the database of studies
def uorf_reads(database, transcripts):
    os.chdir('/scratch/09369/umapaul/data') #EDIT

    genes = pd.read_csv(transcripts, sep=',', header=None)[0].values.tolist()
    starts = pd.read_csv(transcripts, sep=',', header=None)[1].values.tolist()
    stops = pd.read_csv(transcripts, sep=',', header=None)[2].values.tolist()

    study_dict, type_dict = study_folder(database)
    gene_reads = {}
    offset_dic = []
    total = len(study_dict)
    current = 1

    for study in study_dict.keys():
        study_path = ("/scratch/09369/umapaul/data/%s" % study) #EDIT
        print(str(current) + "/" + str(total) + ": " + str(study))
        current += 1

        if os.path.exists(study_path):
            os.chdir(study_path)
            for experiment in study_dict[study]:
                experiment_path = (study_path + "/ribo/experiments/%s.ribo" % experiment) # EDIT
                if os.path.exists(experiment_path):
                    r_file = ribo_data(experiment_path)

                    for j in r_file.experiments:  # run through each experiment
                        mmin,mmax,read_pct = intevl(experiment_path, j)
                        offset = psite_offset(r_file, j, mmin, mmax)
                        exp_dic = {}

                        if read_pct < 0.85:
                            continue

                        # df = r_file.get_coverage(experiment=j, range_lower=mmin, range_upper=mmax, alias=True)
                        CDS_reads = get_region_reads(r_file, "CDS", mmin, mmax, j).add(
                            get_region_reads(r_file, "UTR3_junction", mmin, mmax, j).add(
                                get_region_reads(r_file, "UTR5_junction", mmin, mmax, j)))

                        exp_reads = {}
                        frame1_reads = {}
                        frame2_reads = {}
                        frame3_reads = {}


                        for k in range(mmin, (mmax + 1)):

                            df = r_file.get_coverage(experiment=experiment, range_lower=k, range_upper=k,
                                                          alias=True)

                            add_dic(j, (k, offset[k]), exp_dic)

                            for l in range(len(genes)):  # run through the gene names
                                gene = genes[l]
                                start = starts[l]
                                stop = stops[l]

                                if (start <= offset[k]):
                                    continue

                                # transcript_info.write(str(ribo_object.transcript_lengths))
                                try:
                                    coverage = df[gene]  # store the coverage info of the gene
                                    offset_reads = np.sum(coverage[start - offset[k]: stop - offset[k]])


                                    if (gene, start, stop) in exp_reads.keys():
                                        exp_reads[(gene, start, stop)] += offset_reads
                                        # inframe_reads[(gene, start, stop)] += inframe_read
                                        frame1_reads[(gene, start, stop)] += np.sum(coverage[start - offset[k]: stop - offset[k]: 3])
                                        frame2_reads[(gene, start, stop)] += np.sum(coverage[start - offset[k] + 1: stop - offset[k]: 3])
                                        frame3_reads[(gene, start, stop)] += np.sum(coverage[start - offset[k] + 2: stop - offset[k]: 3])
                                    else:
                                        exp_reads[(gene, start, stop)] = offset_reads
                                        frame1_reads[(gene, start, stop)] = np.sum(
                                            coverage[start - offset[k]: stop - offset[k]: 3])
                                        frame2_reads[(gene, start, stop)] = np.sum(
                                            coverage[start - offset[k] + 1: stop - offset[k]: 3])
                                        frame3_reads[(gene, start, stop)] = np.sum(
                                            coverage[start - offset[k] + 2: stop - offset[k]: 3])

                                    if (k == mmax):
                                        add_dic('Study', study, gene_reads)
                                        add_dic('Experiment', j, gene_reads)
                                        add_dic('Cell_Line', type_dict[j], gene_reads)
                                        add_dic('Gene', gene, gene_reads)
                                        add_dic('uORF Reads', exp_reads[(gene, start, stop)], gene_reads)
                                        add_dic('CDS Reads', CDS_reads.loc[gene][j], gene_reads)
                                        add_dic('Experiment Reads', r_file.info['experiments'][j]['Reads'], gene_reads)
                                        add_dic('Start', start, gene_reads)
                                        add_dic('Stop', stop, gene_reads)
                                        add_dic('Frame1', frame1_reads[(gene, start, stop)],
                                                 gene_reads)
                                        add_dic('Frame2', frame2_reads[(gene, start, stop)],
                                                gene_reads)
                                        add_dic('Frame3', frame3_reads[(gene, start, stop)],
                                                gene_reads)

                                        # add_dic('In_read reads',
                                            #         np.sum(gene_cov[starts[i] - offset[k]: stops[i] - offset[k]: 3]),
                                            #         gene_reads)
                                            # add_dic('Total reads',
                                            #         np.sum(gene_cov[starts[i] - offset[k]: stops[i] - offset[k]]),
                                            #         gene_reads

                                except(KeyError):
                                    continue

                        offset_dic.append(exp_dic)

    return gene_reads, offset_dic

####writes the analysis results to a CSV file
# inputs:
#   database:
def write_output(database='/scratch/09369/umapaul/data/mouse_filtered_complete.csv',
                 transcripts='/scratch/09369/umapaul/data/ltdstart_gene_codons.csv',
                 outfile='/scratch/09369/umapaul/output/ltdstart_output_ribotish.csv',
                 outfile_offsets = "/scratch/09369/umapaul/output/offset_ltdstart.csv"):
    output, offset = uorf_reads(database, transcripts)
    # to write values as columns
    output = pd.DataFrame.from_dict(output)
    output.to_csv(outfile, index=False)

    with open(outfile_offsets, "w", newline="") as file:
        writer = csv.writer(file)

        # Write the header row
        writer.writerow(["Experiment", "Read Length", "P-site Offset"])
        for exp in offset:
            # Write the data rows
            for experiment, data in exp.items():
                for read_length, p_site_offset in data:
                    writer.writerow([experiment, read_length, p_site_offset])

write_output()
#study_folder('/scratch/09369/umapaul/data/fin_version_database_mouse_only.csv')


