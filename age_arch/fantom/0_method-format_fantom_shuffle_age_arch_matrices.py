#!/usr/bin/env python
# coding: utf-8

# 20200924
# sarahfong
# after running the age_enhancers pipeline
# format FANTOM enhancer and 100x matched shuffle datasets
# for analysis of age architecture features

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


path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/"
samples = glob.glob("%sall_unique_fantom_erna_112_tissue.bed" % path)


shuffle_path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/shuffle/breaks/"
shuffle_files = glob.glob("%sshuf-all_fantom_enh_112_tissues-*_age_breaks.bed" % shuffle_path)
shuffle_files

# # load the genomic background

#%% In[4]:


syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2", "mya", "mya2"]] # whittle down the df


#%% cell_line/ tissue descriptions


desc_file = "/dors/capra_lab/data/fantom/fantom5/facet_expressed_enhancers/sample_id_descriptions.txt"
desc_df= pd.read_csv(desc_file, sep = '\t', header = None)

# #%% Initial load and merge all FANTOM eRNA x TFBS datasets together

#%% In[11]: Function for formatting fantom/shuffle dataframes


def formatDF(df, datatype):


    df = df.loc[df.chr_syn !="chrX"] # autosomes only

    df["mrca"] = df["mrca"].round(3) # round mrca value
    df["enh_len"]= df.end_enh - df.start_enh

    if "chr_syn" in list(df):

        df["syn_id"] = df["chr_syn"] + ":" + df["start_syn"].map(str) + "-" + df["end_syn"].map(str)
        df["syn_len"]= df.end_syn - df.start_syn
        df["seg_rep"] = df["syn_len"].divide(df["enh_len"]) # percent of the enhancer occupied by each syntenic block

        #####FILTER out short syntenic blocks#####

        df = df.loc[df["syn_len"]>=6]

        #####ARCHITECTURE#####

        df["code"] = "complex_core"
        df.loc[(df["core"] ==0),"code"] = "derived"
        df.loc[(df["core_remodeling"] ==0), "code"] = "simple"

        df["seg_index"] =  (df["seg_index"] +1) # pseudocount for measuring break density

    df["arch"] = 0
    df.loc[(df["core_remodeling"] ==0),"arch"] = "simple"
    df.loc[(df["core_remodeling"] ==1), "arch"] = "complexenh"

    df = pd.merge(df, syn_gen_bkgd,\
     how = "left", on = "mrca" ).sort_values(by="mrca_2")# add species annotation
    df["mrca_2"] = df["mrca_2"].round(3) # round mrca value

    df = df.drop_duplicates()

    # get the max breaks and oldest mrca
    breakDF = df.groupby(["enh_id", "chr_enh","start_enh", "end_enh",\
     "core_remodeling", "arch"])[["seg_index", "mrca", "enh_len"]].max().reset_index()
    breakDF["mrca"] = breakDF["mrca"].round(3) # round mrca value
    breakDF = pd.merge(breakDF, syn_gen_bkgd,\
    how = "left", on = "mrca").sort_values(by="mrca_2")# add species annotation

    breakDF["seg_den"] = (breakDF["seg_index"]).divide(breakDF["enh_len"])

    breakDF["datatype"] = datatype

    breakDF = breakDF.drop_duplicates()

    return df, breakDF

#%% In[15]:# # format 100 shuffles dataframes

shuf_dict={}
shuf_break_dict={}
for f in shuffle_files:

    shuffle = pd.read_csv(f, sep = '\t', header = None, low_memory = False)
    if "oldest_age" in f:
        shuffle.columns = ["chr_enh", "start_enh", "end_enh",\
        "enh_id","mrca", "enh_len", "seg_index", "core_remodeling"]# rename columns
    else:
        shuffle.columns = ["chr_syn", "start_syn","end_syn","enh_id",
                       "chr_enh", "start_enh", "end_enh", "seg_index",
                       "core_remodeling", "core", "mrca"]
    shuf_id = (f.split("/")[-1]).split(".")[0]

    formatted_shuffle, formatted_shuffle_breaks = formatDF(shuffle, "shuffle") # format df
    formatted_shuffle["shuf_id"] = shuf_id
    shuf_dict[shuf_id] = formatted_shuffle # add to dict

    shuf_break_dict[shuf_id] = formatted_shuffle_breaks # add to dict

    shuffle = formatted_shuffle.loc[formatted_shuffle.enh_len<10000]


     # # Write summary shuffle enhancer file w/ oldest age/ most breaks/arch

        # Don't remove the longest enhancer in the shuffle dataset. This matches FANTOM.
    shufBreaks = formatted_shuffle_breaks.loc[formatted_shuffle_breaks.enh_len<10000]
    shufBreaks = shufBreaks[['chr_enh', 'start_enh', 'end_enh', 'enh_id', 'core_remodeling',
     'arch', 'seg_index',
     'mrca', 'enh_len', 'taxon','mrca_2', 'taxon2', 'mya', 'mya2', 'seg_den', 'datatype']]

    out_enh_shufBreaks = "%sSHUFFLE_FANTOM_%s_summary_matrix.bed" % (shuffle_path, shuf_id)

    shufBreaks.to_csv(out_enh_shufBreaks, sep = '\t', index = False, header = False)

shuffle = pd.concat(shuf_dict.values())
shuffle.head()

#%% In[17]:# # FANTOM enhancer architecture
# open combined eRNA file of 112 FANTOM samples
# sample_id assignment
# sample_desc assignment
# add dataframe to dictionary
# concatenate dataframes

enh_dict = {} # {sample_id: enh_tfbs_density df}
enh_break_dict = {}
for sample in samples:
    if os.path.getsize(sample)>0:
        df = pd.read_csv(sample, sep='\t', header = None)
        df = df[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]].drop_duplicates()

        # rename the columns
        df.columns = ["chr_syn", "start_syn","end_syn","enh_id",
                       "chr_enh", "start_enh", "end_enh", "seg_index",
                       "core_remodeling", "core", "mrca"]

        sample_id = "".join((sample.split("/")[-1]).split("_")[0])

        df["id"] = sample_id # add sample_id to column

        formatted_df, formatted_df_breaks = formatDF(df, "FANTOM") # format df

        enh_dict[sample_id] = formatted_df # add sample df to dictionary

        enh_break_dict[sample_id] = formatted_df_breaks

df = pd.concat(enh_dict.values())

df.head()

#%% In[19]: write the full FANTOM enhancer matrix.

# ## remove large random FANTOM enhancer. This must be an error.
print(df.enh_len.max())

final_merge = df.loc[df.enh_len !=df.enh_len.max()] # remove this weirdo

print(final_merge.enh_len.max())

# save the file
out_enh = "%sFANTOM_enh_age_arch_full_matrix.tsv" % path
final_merge.to_csv(out_enh, sep = '\t', index = False)

final_merge.head()

#%% In[33]:
# # Write summary FANTOM enhancer file w/ oldest age/ most breaks/arch

breaks = pd.concat(enh_break_dict.values())
print(breaks.shape)

# remove the longest enhancer in the FANTOM dataset...

breaks = breaks.loc[breaks.enh_len <10000] # remove this weirdo

out_enh_breaks = "%sFANTOM_enh_age_arch_summary_matrix.tsv" % path
breaks["shuf_id"] = "FANTOM"

breaks = breaks[['chr_enh', 'start_enh', 'end_enh', 'enh_id', 'core_remodeling',
 'arch', 'seg_index',
 'mrca', 'enh_len', 'taxon','mrca_2', 'taxon2', 'mya', 'mya2', 'seg_den', 'datatype']]

breaks.to_csv(out_enh_breaks, sep = '\t', index = False)

# make a bedfile
out_enh_breaks_bed = "%sFANTOM_enh_age_arch_summary_matrix.bed" % path
breaks.to_csv(out_enh_breaks_bed, sep = '\t', index = False, header = False)
breaks.head()
#%%
list(breaks)
#%% In[32]:# # Write shuffle full dataframe

shuffle = shuffle.loc[shuffle.enh_len<10000]
print(shuffle.shape)
out_shuf = "%sSHUFFLED_FANTOM_enh_age_arch_full_matrix.tsv" % path
shuffle.to_csv(out_shuf, sep = '\t', index = False)
shuffle.head()

#%% # # Write summary shuffle enhancer file w/ oldest age/ most breaks/arch

shufBreaks = pd.concat(shuf_break_dict.values())
print(shufBreaks.shape)

# Don't remove the longest enhancer in the shuffle dataset. This matches FANTOM.
shufBreaks = shufBreaks.loc[shufBreaks.enh_len<10000]
shufBreaks = shufBreaks[['chr_enh', 'start_enh', 'end_enh', 'enh_id', 'core_remodeling',
 'arch', 'seg_index',
 'mrca', 'enh_len', 'taxon','mrca_2', 'taxon2', 'mya', 'mya2', 'seg_den', 'datatype']]

out_enh_shufBreaks = "%sSHUFFLE_FANTOM_enh_age_arch_summary_matrix.tsv" % path

shufBreaks.to_csv(out_enh_shufBreaks, sep = '\t', index = False)

shufBreaks.head()

#%% make bed file of shuffled.
out_shuf_breaks_bed = "%sSHUFFLE_FANTOM_enh_age_arch_summary_matrix.bed" % path
shufBreaks.to_csv(out_shuf_breaks_bed, sep = '\t', index = False, header = False)
shufBreaks.head()
