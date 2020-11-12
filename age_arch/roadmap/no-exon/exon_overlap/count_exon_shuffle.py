
import glob
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

RE = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/results/for_publication/non-genic/"
path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/shuffle/breaks/"



noexon_fs = glob.glob("%sno-exon_E*.bed" % path)
wexon_fs = glob.glob("%sexonOverlap_E*.bed" % path)


noexon_fs = glob.glob("%sno-exon_E*.bed" % path)
wexon_fs = glob.glob("%sexonOverlap_E*.bed" % path)


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
#%%
desc_path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/"
desc_df = pd.read_csv("%sroadmap_hg19_sample_id_desc.csv"% desc_path, header = None)
desc_df.columns = ["sid", "desc"]
desc_df.head()
#%%



sid_dict = {}



for f in noexon_fs:

    sid = (f.split("/")[-1]).split("_")[1]


    if sid not in sid_dict:
        print(sid)

        new_len = file_len(f)

        sid_dict[sid] = [new_len] # add


#%%

for f in wexon_fs:

    sid = (f.split("/")[-1]).split("_")[1]

    if sid in sid_dict.keys():

        new_list = sid_dict[sid] # get the list

        new_len = file_len(f)


        new_list.append(new_len)


        sid_dict[sid] = new_list # add
#%%



val = 0
df = pd.DataFrame()
for key, val in sid_dict.items():
    newdf = pd.DataFrame({"sid":[key],
    "no_ex_enh":[val[0]],
    'ex_enh':[val[1]],
    'total_enh': [(val[0]+val[1])]})
    newdf["frac_noex"] = newdf.no_ex_enh.divide(newdf.total_enh)
    newdf["frac_ex"] = newdf.ex_enh.divide(newdf.total_enh)
    if val == 0:
        df = newdf.copy()
    else:
        df = pd.concat([df, newdf])
    val=+1

#%%
df.head()
df.describe()

df.sort_values(by = "frac_ex", ascending = False)


# add tissue description
df.sid = df.sid.apply(lambda x: x.split(".")[0])
df = pd.merge(df, desc_df, how = "left", on  ="sid")
df["sid2"] = df.sid +"-"+df.desc



#%%
df.head()
#%%
df.describe()
#%%

fig, ax = plt.subplots(figsize = (6,10))

sns.barplot(y= "sid", x = "frac_ex",
data = df.sort_values(by = "frac_ex", ascending = False))
#ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
plt.savefig("%sexon_overlap_fraction_shuffle.pdf" %RE, bbox_inches = "tight")

#%%
fig, ax = plt.subplots(figsize = (6,10))

sns.barplot(y= "sid", x = "ex_enh",
data = df.sort_values(by = "frac_ex", ascending = False))
#ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
plt.savefig("%sexon_overlap_count_shuffle.pdf" %RE, bbox_inches = "tight")
#%%
df.head()
df.sid2 = df.sid.apply(lambda x: x.split("-")[0])
#%%
df.sid2.unique()
