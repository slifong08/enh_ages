import glob
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/breaks/"

noexon_fs = glob.glob("%snon-genic/no-exon_ROADMAP_E*_enh_and_shuf_age_arch_summary_matrix.bed" % path)
wexon_fs = glob.glob("%sROADMAP_E*_enh_and_shuf_age_arch_summary_matrix.bed" % path)


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


#%%



sid_dict = {}



for f in noexon_fs:

    sid = (f.split("/")[-1]).split("_")[2]

    if sid not in sid_dict:

        df = pd.read_csv(f, sep ='\t', header = None, usecols = [3,16])
        shuf = df.loc[df[16].str.contains("shuf")]
        enh = df.loc[~df[16].str.contains("shuf")]

        new_len = len(enh)
        new_shuf = len(shuf)



        sid_dict[sid] = [new_len, new_shuf] # add

#%%
sid_dict
#%%

for f in wexon_fs:

    sid = (f.split("/")[-1]).split("_")[1]

    if sid in sid_dict.keys():

        new_list = sid_dict[sid] # get the list

        df = pd.read_csv(f, sep ='\t', header = None, usecols = [3,16], error_bad_lines = False)
        shuf = df.loc[df[16].str.contains("shuf")]
        enh = df.loc[~df[16].str.contains("shuf")]

        new_len = len(enh)
        new_shuf = len(shuf)
        new_list.append(new_len)
        new_list.append(new_shuf)

        sid_dict[sid] = new_list # add
#%%
len(df[16].unique())
sid
#%%
val = 0
df = pd.DataFrame()
for key, val in sid_dict.items():
    newdf = pd.DataFrame({"sid":[key],
    "no_ex_enh":[val[0]],
    "no_ex_shuf":[val[1]],
    'total_enh':[val[2]],
    'total_enh':[val[3]],
    })
    newdf["delta"] = newdf.total_enh - newdf.no_ex_enh
    if val == 0:
        p
        df = newdf.copy()
    else:
        df = pd.concat([df, newdf])
    val=+1
#%%
df.head()
#%%
df.sort_values(by = "delta", ascending = False)
df["log2_ex_overlap"] = np.log10(df.delta)
#%%
fig, ax = plt.subplots(figsize = (6,30))

sns.barplot(y= "sid", x = "log2_ex_overlap",
data = df.sort_values(by = "delta", ascending = False))
#ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
