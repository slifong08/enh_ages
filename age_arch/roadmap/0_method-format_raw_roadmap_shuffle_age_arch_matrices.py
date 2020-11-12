#!/usr/bin/env python
# coding: utf-8

# 20200924
# sarahfong
# after running the age_enhancers pipeline
# format RANDOM enhancer and 100x matched shuffle datasets
# for analysis of age architecture features

###########################################################################
#
#########################################################################



#%% In[1]:

import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os, sys
import datetime
last_run = datetime.datetime.now()
today = (datetime.date.today())
print("last run", datetime.datetime.now())

#%% In[3]: load paths


path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/"

samples = glob.glob("%s**/*age_breaks.bed" % path, recursive = True)
samples = glob.glob("%sshuffle/breaks/*age_breaks.bed" % path)
print(len(samples))

already_done = glob.glob("%sbreaks/ROADMAP_*_enh_summary_matrix.bed" % path)
already_done_list = []

for i in already_done:
    sid = (i.split("/")[-1]).split("_")[1]
    already_done_list.append(sid)

print((already_done_list), len(already_done_list))



to_do_list = ["E055", "E111", "E078", "E102", "E103", "E008", "E098", "E012",\
 "E114", "E056", "E105", "E034", "E066", "E029", "E003", "E063"]


#%% # load the genomic background


syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2", "mya", "mya2"]] # whittle down the df


#%% cell_line/ tissue descriptions

## CHANGE PATH ##

desc_file = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/roadmap_hg19_sample_id_desc.csv"
desc_df= pd.read_csv(desc_file, sep = '\t', header = None)

# #%% Initial load and merge all ROADMAP and shuffled (raw and trimmed)

# Function for formatting ROADMAP/shuffle dataframes


def formatDF(df, datatype):


    df = df.loc[df.chr_enh !="chrX"] # autosomes only

    df[[ "start_enh", "end_enh", "seg_index", "core_remodeling"]]\
     = df[[ "start_enh", "end_enh", "seg_index", "core_remodeling"]].astype(int)

    df["mrca"] = df["mrca"].astype(float).round(3) # round mrca value

    df["enh_len"]= df.end_enh - df.start_enh
    df = df.loc[df.enh_len <10000]

    if "chr_syn" in list(df):
        df[["start_syn", "end_syn"]] =  df[["start_syn", "end_syn"]].astype(int)

        df["syn_id"] = df["chr_syn"] + ":" + df["start_syn"].map(str) + "-" + df["end_syn"].map(str)
        df["syn_len"]= df.end_syn - df.start_syn
        df["seg_rep"] = df["syn_len"].divide(df["enh_len"]) # percent of the enhancer occupied by each syntenic block

        #####FILTER out short syntenic blocks#####

        df = df.loc[df["syn_len"]>=6] # this was done when running the breaks, but just to make sure.

        #####ARCHITECTURE#####

        df["code"] = "complex_core"
        df.loc[(df["core"] ==0),"code"] = "derived"
        df.loc[(df["core_remodeling"] ==0), "code"] = "simple"

        df["seg_index"] =  (df["seg_index"] +1) # pseudocount for measuring break density

    df["arch"] = 0
    df.loc[(df["core_remodeling"] ==0),"arch"] = "simple"
    df.loc[(df["core_remodeling"] ==1), "arch"] = "complexenh"

    # add grouped ages

    df = pd.merge(df, syn_gen_bkgd,\
     how = "left", on = "mrca" ).sort_values(by="mrca_2")# add species annotation
    df["mrca_2"] = df["mrca_2"].round(3) # round mrca value

    df = df.drop_duplicates()


    # get the max breaks and oldest mrca

    breakDF = df.groupby(["enh_id", "chr_enh","start_enh", "end_enh", "shuf_id", \
     "core_remodeling", "arch"])[["seg_index", "mrca", "enh_len"]].apply(max).reset_index()

    breakDF["mrca"] = breakDF["mrca"].round(3) # round mrca value

    breakDF = pd.merge(breakDF, syn_gen_bkgd,\
    how = "left", on = "mrca").sort_values(by="mrca_2")# add species annotation

    breakDF["seg_den"] = (breakDF["seg_index"]).divide(breakDF["enh_len"])

    breakDF["datatype"] = datatype

    breakDF = breakDF.drop_duplicates()

    return df, breakDF


def getSid(sample_id):
    if "shuf-E" in sample_id: #"shuf-E" in sample_id
        base = sample_id
        sid = sample_id.split("_")[0]

    elif "shuf-trimmed-310_" in sample_id:
        base = sample_id
        sid = sample_id.split("_")[1]


    elif "trimmed-310-E" in sample_id: #
        base = sample_id
        sid =  (base.split("_")[0]).split("-")[2]


    else:
        base = sample_id
        sid = base.split("_")[0]

    return sid


#%% In[17]:# # ROADMAP enhancer architecture
# open ChIP-seq file of 98 ROADMAP samples
# sample_id assignment
# sample_desc assignment
# add dataframe to dictionary
# concatenate dataframes


len(samples)
#%%

raw_sample = {} # make a dict of sid : [sample_file1, shufsample_file1, shufsamplefile2, etc.]
trimmed_sample = {}

for sample in samples:


    if os.path.getsize(sample)>0:
        sample_id = "".join((sample.split("/")[-1]).split(".")[0])
        sid = getSid(sample_id)

        if "trimmed" in sample:
            if sid not in trimmed_sample.keys():
                trimmed_sample[sid] = [sample]
            else:
                sample_list = trimmed_sample[sid]
                sample_list.append(sample)
        else:
            if sid not in raw_sample.keys():
                raw_sample[sid] = [sample]
            else:
                sample_list = raw_sample[sid]
                sample_list.append(sample)
#%%
raw_sample.values()
#%%
for sid, sample_list in raw_sample.items():


    #if sid not in already_done_list:
    print("working on", sid)
    enh_dict = {} # full df
    enh_break_dict = {} # enh summar df

    for sample in sample_list:

        sample_id = "".join((sample.split("/")[-1]).split(".")[0])
        print(sample_id)

        cols = pd.read_csv(sample, sep='\t', nrows = 1).columns

        if "start_enh" in cols:
            df = pd.read_csv(sample, sep='\t').drop_duplicates()
        else:
            df = pd.read_csv(sample, sep='\t', header = None).drop_duplicates()

        num_cols = len(list(df))
        print(num_cols)
        print(df.head())
        if num_cols ==11:
            df.columns = ["chr_syn", "start_syn","end_syn","enh_id",
                               "chr_enh", "start_enh", "end_enh", "seg_index",
                               "core_remodeling", "core", "mrca"]
            df["shuf_id"] = sample_id
        if num_cols == 9:
            df.columns = ["chr_enh", "start_enh", "end_enh", "old_len", "enh_id",
        	"seg_index", "core_remodeling", "arch", "mrca"]
            df = df.loc[df.mrca != "max_age"] # these are header columns that need to get removed.

        elif "shuf-trimmed" in sample_id:
            df.columns = ["chr_syn", "start_syn","end_syn", "shuf_id", "enh_id",
                       "chr_enh", "start_enh", "end_enh", "seg_index",
                       "core_remodeling", "core", "mrca"]
        elif "sid" in list(df):
            df.columns = ["chr_syn", "start_syn","end_syn","enh_id",
                           "chr_enh", "start_enh", "end_enh", "seg_index",
                           "core_remodeling", "core", "mrca", "shuf_id", "enh_len"]
        elif "chr_syn" in cols and num_cols !=11:

            df.columns = ["chr_syn", "start_syn","end_syn","enh_id",
                               "chr_enh", "start_enh", "end_enh", "seg_index",
                               "core_remodeling", "core", "mrca", "enh_len",  "shuf_id"]
        #else:
        #    df.columns = ["chr_syn", "start_syn","end_syn","enh_id",
        #                   "chr_enh", "start_enh", "end_enh", "seg_index",
        #                   "core_remodeling", "core", "mrca",  "enh_len", "shuf_id",]

        df["id"] = sample_id # add sample_id to column

        ## CHANGE sid in format ##
        formatted_df, formatted_df_breaks = formatDF(df, sample_id) # format df

        enh_dict[sample_id] = formatted_df # add sample df to dictionary

        enh_break_dict[sample_id] = formatted_df_breaks

    df = pd.concat(enh_dict.values()) # concat enhancer and shuffles together

    print(df.shape)

    final_merge = df.loc[df.enh_len<10000]

    print(final_merge.shape)

    print(final_merge.enh_len.max())

    # save the file
    out_enh = "%sbreaks/ROADMAP_%s_enh_and_shuf_age_arch_full_matrix.tsv" % (path, sid)
    #final_merge.to_csv(out_enh, sep = '\t', index = False)

    # # Write summary ROADMAP enhancer file w/ oldest age/ most breaks/arch

    breaksAll = pd.concat(enh_break_dict.values())

    breaks = breaksAll.loc[breaksAll.enh_len<10000]

    out_enh_breaks = "%sshuffle/breaks/ROADMAP_%s_enh_summary_matrix.tsv" % (path, sid)

    breaks = breaks[['chr_enh', 'start_enh', 'end_enh', "shuf_id", 'enh_id', 'core_remodeling',
    'arch', 'seg_index',
    'mrca', 'enh_len', 'taxon','mrca_2', 'taxon2', 'mya', 'mya2', 'seg_den', 'datatype']]

    #breaks.to_csv(out_enh_breaks, sep = '\t', index = False)

    # make a bedfile
    out_enh_breaks_bed = "%sbreaks/ROADMAP_%s_enh_summary_matrix.bed" % (path, sid)
    breaks.to_csv(out_enh_breaks_bed, sep = '\t', index = False, header = False)

#%%
print(sample_id)
#%%
df = pd.read_csv(sample, sep='\t', header = None).drop_duplicates()
df.head()
#%%
cols = pd.read_csv(sample, sep='\t', nrows = 1).columns
len(cols)
