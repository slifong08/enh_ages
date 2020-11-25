#!/usr/bin/env python
# coding: utf-8

# 20200924
# sarahfong
# after running the age_enhancers pipeline
# format RANDOM enhancer and 100x matched shuffle datasets
# for analysis of age architecture features

###########################################################################
# STOP
# CHANGE ROADMAP TO ROADMAP PATHS BEFORE running
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

trimmed_path = "%strimmed/breaks/" % path
## CHANGE PATH ##
samples = glob.glob("%s**/*parallel_breaks.bed" % path, recursive = True)
samples = glob.glob("%strimmed-310-E*_age_breaks.bed" % trimmed_path, recursive = True)


# # load the genomic background

#%% In[4]:


syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2", "mya", "mya2"]] # whittle down the df


#%% cell_line/ tissue descriptions

## CHANGE PATH ##

desc_file = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/roadmap_hg19_sample_id_desc.csv"
desc_df= pd.read_csv(desc_file, sep = '\t', header = None)

# #%% Initial load and merge all ROADMAP eRNA x TFBS datasets together

#%% In[11]: Function for formatting ROADMAP/shuffle dataframes


def formatDF(df, datatype):


    df = df.loc[df.chr_enh !="chrX"] # autosomes only

    df[[ "start_enh", "end_enh", "seg_index", "core_remodeling"]]\
     = df[[ "start_enh", "end_enh", "seg_index", "core_remodeling"]].astype(int)

    df["mrca"] = df["mrca"].astype(float).round(3) # round mrca value

    df["enh_len"]= df.end_enh - df.start_enh
    df = df.loc[df.enh_len <10000]

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
    if "old_len" in list(df):
        breakDF = df.groupby(["enh_id", "chr_enh","start_enh", "end_enh", "old_len", \
         "core_remodeling", "arch"])[["seg_index", "mrca", "enh_len"]].max().reset_index()
    else:
        breakDF = df.groupby(["enh_id", "chr_enh","start_enh", "end_enh", \
     "core_remodeling", "arch"])[["seg_index", "mrca", "enh_len"]].max().reset_index()

    breakDF["mrca"] = breakDF["mrca"].round(3) # round mrca value

    breakDF = pd.merge(breakDF, syn_gen_bkgd,\
    how = "left", on = "mrca").sort_values(by="mrca_2")# add species annotation

    breakDF["seg_den"] = (breakDF["seg_index"]).divide(breakDF["enh_len"])

    breakDF["datatype"] = datatype

    breakDF = breakDF.drop_duplicates()

    return df, breakDF


#%% In[17]:# # ROADMAP enhancer architecture
# open ChIP-seq file of 98 ROADMAP samples
# sample_id assignment
# sample_desc assignment
# add dataframe to dictionary
# concatenate dataframes

enh_dict = {} # {sample_id: enh_tfbs_density df}
enh_break_dict = {}

for sample in samples:
    sid = (sample.split("/")[-1]).split("_")[0]

    if os.path.getsize(sample)>0:
        df = pd.read_csv(sample, sep='\t')

        num_cols = len(list(df))

        if num_cols==11:
            # rename the columns
            df.columns = ["chr_syn", "start_syn","end_syn","enh_id",
                       "chr_enh", "start_enh", "end_enh", "seg_index",
                       "core_remodeling", "core", "mrca"]

        elif num_cols == 9:
            df.columns = ["chr_enh", "start_enh", "end_enh", "old_len", "enh_id",
            	"seg_index", "core_remodeling", "arch", "mrca"]

        df = df.loc[df.mrca != "max_age"] # these are header columns that need to get removed.

        sample_id = "".join((sample.split("/")[-1]).split(".")[0])

        df["id"] = sample_id # add sample_id to column

        ## CHANGE sid in format ##
        formatted_df, formatted_df_breaks = formatDF(df, sample_id) # format df



        breaks = formatted_df_breaks[['chr_enh', 'start_enh', 'end_enh', 'enh_id', 'core_remodeling',
         'arch', 'seg_index',
         'mrca', 'enh_len', 'taxon','mrca_2', 'taxon2', 'mya', 'mya2', 'seg_den', 'datatype']]

        out_enh_breaks_bed = "%s%s_summary_matrix.bed" % (trimmed_path, sid)

        breaks.to_csv(out_enh_breaks_bed, sep = '\t', index = False, header = False)

        enh_dict[sample_id] = formatted_df # add sample df to dictionary

        enh_break_dict[sample_id] = formatted_df_breaks

#%%

df = pd.concat(enh_dict.values())

df.head()
#%%

# ## remove large random ROADMAP enhancer. This must be an error.
print(df.enh_len.max())
print(df.shape)

#%%
final_merge = df.loc[df.id.isin(ids)] # remove this weirdo
final_merge = final_merge.loc[final_merge.enh_len<10000]
print(final_merge.shape)
#%%
print(final_merge.enh_len.max())

# save the file
out_enh = "%sROADMAP_enh_age_arch_full_matrix.tsv" % path
final_merge.to_csv(out_enh, sep = '\t', index = False)

final_merge.head()

#%% In[33]:
# # Write summary ROADMAP enhancer file w/ oldest age/ most breaks/arch

breaksAll = pd.concat(enh_break_dict.values())

#%%

breaks = breaksAll.loc[breaksAll.datatype.isin(ids)] # remove this weirdo
breaks = breaks.loc[breaks.enh_len<10000]

breaks.head()
#%%


out_enh_breaks = "%sROADMAP_enh_age_arch_summary_matrix.tsv" % path

breaks = breaks[['chr_enh', 'start_enh', 'end_enh', 'enh_id', 'core_remodeling',
 'arch', 'seg_index',
 'mrca', 'enh_len', 'taxon','mrca_2', 'taxon2', 'mya', 'mya2', 'seg_den', 'datatype']]

breaks.to_csv(out_enh_breaks, sep = '\t', index = False)

# make a bedfile
out_enh_breaks_bed = "%sROADMAP_enh_age_arch_summary_matrix.bed" % path
breaks.to_csv(out_enh_breaks_bed, sep = '\t', index = False, header = False)
breaks.head()
#%%



#%% In[32]:# # Write shuffle full dataframe
shuffle = df.loc[df.id.str.contains("shuf")]
shuffle = shuffle.loc[shuffle.enh_len<10000]
print(shuffle.shape)
out_shuf = "%sSHUFFLED_ROADMAP_enh_age_arch_full_matrix.tsv" % path
shuffle.to_csv(out_shuf, sep = '\t', index = False)
shuffle.head()

#%% # # Write summary shuffle enhancer file w/ oldest age/ most breaks/arch

shufBreaks = breaksAll.loc[breaksAll.datatype.str.contains("shuf")]
shufBreaks = shufBreaks.loc[shufBreaks.enh_len<10000]
print(shufBreaks.shape)

# Don't remove the longest enhancer in the shuffle dataset. This matches ROADMAP.
print(shufBreaks.enh_len.max(), "is this a long weirdo eRNA?")

# rearrange the columns as if it were a .bed file
shufBreaks = shufBreaks[['chr_enh', 'start_enh', 'end_enh', 'enh_id', 'core_remodeling',
 'arch', 'seg_index',
 'mrca', 'enh_len', 'taxon','mrca_2', 'taxon2', 'mya', 'mya2', 'seg_den', 'datatype']]
out_enh_shufBreaks = "%sSHUFFLE_ROADMAP_enh_age_arch_summary_matrix.tsv" % path

shufBreaks.to_csv(out_enh_shufBreaks, sep = '\t', index = False)
# make a bedfile
out_enh_shufBreaks_bed = "%sSHUFFLE_ROADMAP_enh_age_arch_summary_matrix.bed" % path
shufBreaks.to_csv(out_enh_shufBreaks_bed, sep = '\t', index = False, header = False)

shufBreaks.head()


#%% individual dataset export

for key in ids:

    print(key)
    final_merge = df.loc[df.id == key] # remove this weirdo
    final_merge = final_merge.loc[final_merge.enh_len<10000]
    print(final_merge.shape)

    print(final_merge.enh_len.max())

    # save the file
    out_enh = "%sROADMAP_%s_enh_age_arch_full_matrix.tsv" % (path, key)
    final_merge.to_csv(out_enh, sep = '\t', index = False)

    ### SUMMARY MATRIX
    breaks = breaksAll.loc[breaksAll.datatype == key] # remove this weirdo
    breaks = breaks.loc[breaks.enh_len<10000]


    out_enh_breaks = "%sROADMAP_%s_summary_matrix.tsv" % (path, key)

    if "old_len" in list(breaks):
        breaks = breaks[['chr_enh', 'start_enh', 'end_enh', "old_len", 'enh_id', 'core_remodeling',
    'arch', 'seg_index',
    'mrca', 'enh_len', 'taxon','mrca_2', 'taxon2', 'mya', 'mya2', 'seg_den', 'datatype']]
    else:
        breaks = breaks[['chr_enh', 'start_enh', 'end_enh', 'enh_id', 'core_remodeling',
'arch', 'seg_index',
'mrca', 'enh_len', 'taxon','mrca_2', 'taxon2', 'mya', 'mya2', 'seg_den', 'datatype']]
    breaks.to_csv(out_enh_breaks, sep = '\t', index = False)

    # make a bedfile
    out_enh_breaks_bed = "%sROADMAP_%s_enh_age_arch_summary_matrix.bed" % (path, key)
    breaks.to_csv(out_enh_breaks_bed, sep = '\t', index = False, header = False)
    breaks.head()
