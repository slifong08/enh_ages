#!/usr/bin/env python
# coding: utf-8

# In[34]:


# intro

# created 2019-08-05
# sarahfong

# Goal - multiintersect all roadmap histone enhancers with a count for the number of overlapping coordinates between datasets

# Method - multiintersect

## step 1 - make bed file with all non-redundant enhancers, named "all_fantom_enhancers". This will be the -a arguement
## step 2 - bedtools intersect -a [all_roadmap_enhancers] -b [E005_core_ages.bed, E004_core_ages.bed ... etc.]

## -c (count number of overlapping features)

import glob
import os, sys
import pandas as pd
import datetime
import seaborn as sns
import matplotlib.pyplot as plt
#%%


trim_len = "mean"


test_list = [] # collect sids that have been trimmed
mean_len_dict = {} # dict of sid:mean length
relative_simple_dict = {} # dict of sid: relative_simple cutoff
already_done = [] # list of shuffled sids that have already been trimmed and assigned cognate relative simple cutoff

# #  Trim function

# In[14]:


def trim(bedfile, trim_len, relative_simple):

    df = pd.read_csv(bedfile, sep ='\t', header = None, usecols  = [0,1,2,3,5]) # open file

    df.columns = ["chr", "start", "end", "old_enh_id", "seg_index"] # name columns

    #df["old_enh_id"] = df.chr + ":" + df.start.map(str) + "-"+ df.end.map(str)
    df["old_len"]= df.end - df.start # calculate enhancer length

    if relative_simple ==0:
        relative_simple = df.seg_index.median()
        df["core_remodeling"] = 0
        df.loc[df.seg_index >= relative_simple, "core_remodeling"] = 1

    else:
        df["core_remodeling"] = 0
        df.loc[df.seg_index >= relative_simple, "core_remodeling"] = 1

    if trim_len == "mean": # take the mean

        trim_len = df["old_len"].mean().round(0) # mean enhancer length in dataset

    print(df.old_enh_id[0], trim_len)
    df["midpoint"] = (df.start + (trim_len)/2).astype(int) # identify the midpoint of each enhancer

    df["new_len"] = trim_len

    df["start_new"] =((df.midpoint - (trim_len/2)).round(0)).astype(int) # calculate new start as the midpoint - (mean length/2)

    df["end_new"] = ((df.midpoint + (trim_len/2)).round(0)).astype(int)

    trimmed = df[["chr", "start_new", "end_new", "old_enh_id",
    "old_len", "new_len", "seg_index", "core_remodeling"]].drop_duplicates()

    return trimmed, trim_len, relative_simple


# In[ ]:


source_path ="/dors/capra_lab/data/roadmap_epigenomics/release_v9/consolidated/histone/h3k27ac_plus_h3k4me3_minus_peaks/"
glob_path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/Hsap_H3K27ac_plus_H3K4me3_minus_E*/cleaned_Hsap_H3K27ac_plus_H3K4me3_minus_E*/breaks/"
samples = glob.glob("%sno-exon*.bed" % glob_path)

for sample in samples:
    s = ((sample.split("/")[-1]).split("_")[-1]).split(".")[0]
    if s not in test_list:
        test_list.append(s)
        print(s)

        path = "/".join(sample.split("/")[:-1]) +"/"

        outpath = "%strimmed/" %path

        if os.path.exists(outpath) == False:
            os.mkdir(outpath)

        trimmed_df, mean_len, relative_simple = trim(sample, "mean", 0)
        trimmed_df.to_csv("%strimmed-%s-%s.bed" % (outpath, trim_len,  s),
        sep = '\t', header = False, index = False)

        mean_len_dict[s] = mean_len
        relative_simple_dict[s] =relative_simple

#%% Now, do the shuffle based on mean of enhancer dataset
len(mean_len_dict.keys())

done = glob.glob("/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/Hsap_H3K27ac_plus_H3K4me3_minus_E*/cleaned_Hsap_H3K27ac_plus_H3K4me3_minus_E*/shuffle/breaks/trimmed/*.bed")
for d in done:
    sid = ((d.split("/")[-1]).split(".")[1]).split("-")[-1]
    print(sid)
    if sid in mean_len_dict.keys():
        mean_len_dict.pop(sid)
        already_done.append(sid)

#%%
print(len(test_list))
#%%

print(len(mean_len_dict.keys()))
print(len(already_done))
#%%
glob_path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/Hsap_H3K27ac_plus_H3K4me3_minus_E*/cleaned_Hsap_H3K27ac_plus_H3K4me3_minus_E*/shuffle/breaks/"
samples = glob.glob("%sshuf-cleaned_Hsap_H3K27ac_plus_H3K4me3_minus_*_enh_age_arch_summary_matrix.bed" % glob_path)

for sample in samples:
    s = ((sample.split("/")[-4]).split("_")[-1]).split(".")[0]
    path = "/".join(sample.split("/")[:-1]) +"/"
    outpath = "%strimmed/" %path

    if os.path.exists(outpath) == False:
        os.mkdir(outpath)

    if s in mean_len_dict.keys():
        trim_len = mean_len_dict[s]
        relative_simple = relative_simple_dict[s]
        print(s, trim_len, relative_simple)
        trimmed_df, mean_len, relative_simple = trim(sample, trim_len, relative_simple)
        trimmed_df.to_csv("%strimmed-%s-shuf-%s.bed" % (outpath, trim_len,  s),
         sep = '\t', header = False, index = False)
#%%
print(trim_len)

#%%

#%%

#%%
)
