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

fraction_overlap = 0.5
trim_len = "mean"



# #  Trim function

# In[14]:


def trim(bedfile, trim_len):

    df = pd.read_csv(bedfile, sep ='\t', header = None) # open file
    df.columns = ["chr", "start", "end", "id"] # name columns
    df["old_enh_id"] = df.chr + ":" + df.start.map(str) + "-"+ df.end.map(str)
    df["old_len"]= df.end - df.start # calculate enhancer length

    if trim_len == "mean": # take the mean

        trim_len = df["old_len"].mean().round(0) # mean enhancer length in dataset
    print(df.id[0], trim_len)
    df["midpoint"] = (df.start + (trim_len)/2).astype(int) # identify the midpoint of each enhancer

    df["new_len"] = trim_len

    df["start_new"] =((df.midpoint - (trim_len/2)).round(0)).astype(int) # calculate new start as the midpoint - (mean length/2)

    df["end_new"] = ((df.midpoint + (trim_len/2)).round(0)).astype(int)

    trimmed = df[["chr", "start_new", "end_new", "old_enh_id", "old_len", "new_len"]].drop_duplicates()

    return trimmed, trim_len


# In[ ]:


source_path ="/dors/capra_lab/data/roadmap_epigenomics/release_v9/consolidated/histone/h3k27ac_plus_h3k4me3_minus_peaks/"
path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/"
samples = glob.glob("%s*.bed" % source_path)
test_list = []
for sample in samples:
    s = ((sample.split("/")[-1]).split("_")[-1]).split(".")[0]
    test_list.append(s)
    print(s)
# -a

#%%
#test_list =["E050", "E029", "E034", "E069", "E072", "E073", "E118", "E123", "E116"]

outpath = "%strimmed/" %path
if os.path.exists(outpath) == False:
    os.mkdir(outpath)


# In[30]:
sid_mean_len = {}

for test in test_list:

    print(test)
    # Entire Enhancer #
    # Trim enhancers to mean lengths #
    infile = "%sHsap_H3K27ac_plus_H3K4me3_minus_%s.bed" % (source_path, test)

    trimmed_df, mean_len = trim(infile, trim_len)
    trimmed_df.to_csv("%strimmed-%s-%s.bed" % (outpath, trim_len,  test,),
     sep = '\t', header = False, index = False)
    sid_mean_len[test]= mean_len


#%%
df = pd.DataFrame()
for key, val in sid_mean_len.items():
    ndf = pd.DataFrame({"sid":[key], "mean_len":[val]})
    df = pd.concat([df, ndf])

desc_file = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/roadmap_hg19_sample_id_desc.csv"
desc_df= pandas.read_csv(desc_file, header = None)
desc_df.columns = ["sid", "desc"]
df = pd.merge(df, desc_df, how = "left", on ="sid")
df["sid2"] = df.sid + "-" + df.desc
#%%
fig, ax = plt.subplots(figsize = (8,30))
sns.barplot(y = "sid2", x ="mean_len", data = df.sort_values(by = "mean_len"))
RE = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/results/for_publication/"
plt.savefig("%sroadmap_noexon_mean_length.pdf"%RE, bbox_inches = 'tight')
#%%
# In[36]:


source_path = outpath
samples = glob.glob("%strimmed0-*.bed"%source_path)
sample_dict = {}

ITERATIONS = 50
AGE_VAL = 1
BREAK_VAL = 1
TFBS_VAL = 0
SHUF_VAL = 0
RUN_VAL = 1

RUN = 1 # Launch command or dont
SBATCH = 0 # Sbatch or run w/ python interpreter


for sample in samples:
    sample_id = ((sample.split("/")[-1]).split("_")[-1]).split(".")[0]
    sample_dict[sample_id] = sample

outpath
# In[37]:


sample_dict


# In[43]:


print("start", datetime.datetime.now())

for sample_id, file in sample_dict.items():

    if SBATCH ==1:
        cmd = "sbatch /dors/capra_lab/users/fongsl/enh_age/enh_age_git/bin/age_enhancers_w_parallelbreaks.slurm %s %s %s %s %s %s %s" \
                % (file, ITERATIONS, AGE_VAL, BREAK_VAL, TFBS_VAL, SHUF_VAL, RUN_VAL)

        print("SLURM", sample_id, datetime.datetime.now())
        print(cmd)
        if RUN ==1:
            os.system(cmd)
    else:
        cmd = "python /dors/capra_lab/users/fongsl/enh_age/enh_age_git/bin/age_enhancers_w_parallelbreaks.py %s -i %s -a %s -b %s -t %s -sh %s -rt %s"\
                 % (file, ITERATIONS, AGE_VAL, BREAK_VAL, TFBS_VAL, SHUF_VAL, RUN_VAL)

        print("PYTHON", sample_id, datetime.datetime.now())
        print(cmd)
        if RUN ==1:
            os.system(cmd)


print("end", datetime.datetime.now())


# In[42]:


trimmed_df.head()


# In[ ]:
