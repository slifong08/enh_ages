#!/usr/bin/env python
# coding: utf-8

# In[1]:


import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os, sys
from scipy import stats
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')

# sns colors
arch_colors = [ "amber", "dusty purple", "windows blue"]
arch_palette = sns.xkcd_palette(arch_colors)
sns.palplot(arch_palette)

# sns graphing preferences
sns.set(color_codes=True)
sns.set(font_scale=1.5)
sns.set_style("white")
sns.despine(bottom=True, left=True)

import datetime
last_run = datetime.datetime.now()
today = (datetime.date.today())
print("last run", datetime.datetime.now())
RE = "/dors/capra_lab/projects/enhancer_ages/linsight/results/"

AP = RE


# In[2]:


linsight_path = "/dors/capra_lab/projects/enhancer_ages/linsight/data/"
roadmap_fs = glob.glob("%sE*_linsight.bed" % linsight_path)


# In[4]:


syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t')
syn_gen_bkgd["mrca"]=syn_gen_bkgd["mrca"].round(3)
syn_gen_bkgd["mrca_2"]=syn_gen_bkgd["mrca_2"].round(3)
syn_gen_bkgd.head()
syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2"]]

desc_f= "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/roadmap_hg19_sample_id_desc.csv"
desc_df = pd.read_csv(desc_f, header = None)
desc_df.columns = ["sid", "desc"]
desc_df.head()


# In[5]:


def format_df(df, arch_id, desc):
    df.columns = [ 'chr_enh', 'start_enh','end_enh','enh_id',  'core_remodeling', 'seg_index',
                  'mrca', "chr_lin",  "start_lin", "end_lin","linsight_score", "overlap"]
    df['arch'] = arch_id

    df['desc'] = desc

    print(arch_id, desc)

    df['enh_len'] = df.end_enh - df.start_enh


    df.loc[(df.core_remodeling == 0), "arch"] = "simple"
    df.loc[(df.core_remodeling == 1), "arch"] = "complexenh"

    df = df.loc[df.start_lin != -1] # exclude the loci that do not overlap a linsight score

    df.linsight_score = df.linsight_score.astype(float) # turn linsight scores into floats

    df["linsight_id"] = df.chr_lin + ":" + df.start_lin.map(str) +"-"+ df.end_lin.map(str)

    # make a dataframe only based on linsight scores and enhancer architectures

    base_df = df.groupby(["chr_lin", "start_lin", "end_lin", "linsight_score", "arch", "desc"])["mrca"].max().reset_index()

    base_df["lin_len"] = base_df.end_lin - base_df.start_lin

    return base_df
def get_counts(res, strat):

    if strat == 1:
        counts = res.groupby(["mrca_2", "core_remodeling"])["enh_id"].count().reset_index()

        # add empty dataframe for complex human enhancers (do not exist by our def)
        empty = pd.DataFrame({"mrca_2": [0.000],
        "core_remodeling": [1],
        "enh_id": [0]
        })

        counts = pd.concat([counts, empty]) # concat the dataframe

        counts =counts.sort_values(by = ["core_remodeling", 'mrca_2']).reset_index()

    else:
        counts = res.groupby("arch")["enh_id"].count().reset_index()
        counts =counts.sort_values(by = "arch", ascending = False).reset_index()

    counts["enh_id"] = counts["enh_id"].astype(int)
    counts = counts.drop(["index"], axis = 1)

    return counts

# In[6]:

SHUF = 0
linsight_dict = {}
arch_dict = {}
roadmap_f = roadmap_fs[0]

for roadmap_f in roadmap_fs:

    arch_id = (roadmap_f.split("/")[-1]).split("_")[0]# E000 ID
    sid = (roadmap_f.split("/")[-1]).split(".")[0]
    desc = desc_df.desc.loc[desc_df.sid == arch_id].item() # sample description

    val = 0

    sampled_f = f"{linsight_path}/sampled_100k_%s-%s.bed"%(sid, val)
    new_id = arch_id + "-" +str(val)
    print(new_id)

    if SHUF == 1:
        cmd = "shuf -n 10000 %s > %s" % (roadmap_f, sampled_f)
        print(cmd)
        os.system(cmd)
        df = pd.read_csv(sampled_f, sep = '\t', header = None,low_memory=False)

    elif SHUF == 0:
        df = pd.read_csv(sampled_f, sep = '\t', header = None, low_memory=False)



# In[9]:


SHUF = 1
linsight_dict = {}
arch_dict = {}
skip = ["E118", "E050", "E116", "E123"]
for roadmap_f in roadmap_fs:


    arch_id = (roadmap_f.split("/")[-1]).split("_")[0]# E000 ID
    sid = (roadmap_f.split("/")[-1]).split(".")[0]
    desc = desc_df.desc.loc[desc_df.sid == arch_id].item() # sample description

    if arch_id not in skip:
        val = 0
        print(roadmap_f)

        if SHUF == 1:

                df = pd.read_csv(roadmap_f, sep = '\t', header = None,low_memory=False)
                while val<25:
                    new_id = arch_id + "-" +str(val)
                    sampled = df.sample(frac = 0.01)

                    print(arch_id, sid, new_id)
                    base_df = format_df(sampled, arch_id, desc)
                    base_df["sid"]=new_id
                    arch_dict[new_id] = base_df.groupby(["arch", "desc", "sid"])["linsight_score"].describe().reset_index()

                    linsight_dict[new_id] = base_df.groupby(["arch", "mrca", "desc", "sid"])["linsight_score"].describe().reset_index()
                    val +=1
        else:
            fs = glob.glob(f"{linsight_path}sampled_100k_E*_linsight-*.bed")
            for f in fs:
                sid = (f.split("/")[-1]).split(".")[0]
                print(sid)
                df = pd.read_csv(f, sep = '\t', header = None)

                base_df = format_df(df, arch_id, desc)
                base_df["sid"]=sid
                arch_dict[sid] = base_df.groupby(["arch", "desc", "sid"])["linsight_score"].describe().reset_index()
                linsight_dict[sid] = base_df.groupby(["arch", "mrca", "desc", "sid"])["linsight_score"].describe().reset_index()



# In[10]:

arch_dict.values()
linsight_dict.keys()


# In[11]:


linsight = pd.concat(linsight_dict.values())

linsight.mrca = linsight.mrca.round(3)
linsight.head()


# In[12]:


linsight.sid.unique()

linsight.head()
# In[13]:


arch_colors = [ "amber", "faded green"]
arch_palette = sns.xkcd_palette(arch_colors)
sns.palplot(arch_palette)


# In[14]:


fig, ax = plt.subplots(figsize=(8,8))
order = ["simple", "complexenh"]
plot = pd.merge(linsight, syn_gen_bkgd, how = "left")
sns.barplot(x= "taxon2", y = "50%", data =plot.sort_values(by = "mrca_2"),
              hue = "arch", hue_order = order, palette =arch_palette )
ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
ax.set_ylabel("Median LINSIGHT score")
ax.legend(bbox_to_anchor = (1,1))
ax.set_xlabel("")
plt.savefig("%sROADMAP_0.1sampled_median_linsight_per_mrca_no_exon_relative_simple.pdf"%RE, bbox_inches = "tight")

#%%
plot["core_remodeling"] = 0
plot.loc[plot["arch"] == "complexenh", "core_remodeling"] = 1
plot.sort_values(by = ["core_remodeling", "mrca_2"])

counts = plot.groupby(["arch", "mrca_2","core_remodeling"])["count"].mean().reset_index()
empty = pd.DataFrame({"arch":["complexenh"], "mrca_2":[0.000], "core_remodeling": [1], "count":[0]})
counts = pd.concat([counts, empty])
counts = counts.sort_values(by = ["core_remodeling", "mrca_2"]).reset_index()

fig, ax = plt.subplots(figsize=(8,8))
order = ["simple", "complexenh"]
plot = pd.merge(linsight, syn_gen_bkgd, how = "left")
splot = sns.barplot(x= "taxon2", y = "mean", data =plot.sort_values(by = "mrca_2"),
              hue = "arch", hue_order = order, palette =arch_palette )
#STRAT = 0

for n, p in enumerate(splot.patches):
    value = counts.iloc[n]["count"].astype(int)
    splot.annotate(value,
                   (p.get_x() + p.get_width() / 2.,0.01),
                   ha = 'center', va = 'baseline',
                   size=15,
                   rotation = 90,
                   color = "white",
                   xytext = (0, 1),
                   textcoords = 'offset points'
                   )

ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
ax.set_ylabel("Mean LINSIGHT score")
ax.legend(bbox_to_anchor = (1,1))
ax.set_xlabel("")
plt.savefig("%sFIGS32ROADMAP_0.1sampled_mean_linsight_per_mrca_no_exon_relative_simple.pdf"%RE, bbox_inches = "tight")


# In[16]:


RE


# In[17]:


arch = pd.concat(arch_dict.values())
arch["sid2"] = arch.sid.apply(lambda x: "E"+(x.split("E")[1]).split("_")[0])
arch = arch.loc[arch.sid2 != "E050"]
arch = pd.merge(arch, desc_df, how = "left", right_on = "sid", left_on = "sid2")
arch.head()

#%%
simple = arch.loc[arch.arch == "simple", "mean"]
complex = arch.loc[arch.arch == "complexenh", "mean"]

stats.mannwhitneyu(simple, complex)
simple.mean()
complex.mean()
# In[18]:


fig, ax = plt.subplots(figsize=(8,8))
order = ["simple", "complexenh"]

sns.barplot(x= "desc_x", y = "50%", data =arch, hue = "arch",
               hue_order = order, palette =arch_palette,)# showfliers = False)
ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
ax.set_ylabel("Median LINSIGHT score")
ax.set_xlabel("")
ax.legend(bbox_to_anchor = (1,1))
#plt.savefig("%sROADMAP_0.1sampled_median_linsight_overall.pdf"%RE, bbox_inches = "tight")


# In[19]:
arch["core_remodeling"] = 0
arch.loc[arch["arch"] == "complexenh", "core_remodeling"] = 1
arch.sort_values(by = ["core_remodeling"])
counts = arch.groupby(["core_remodeling","arch", "desc_x"])["count"].mean().reset_index()
counts.sort_values(by = "core_remodeling").reset_index()



fig, ax = plt.subplots(figsize=(8,8))
order = ["simple", "complexenh"]

splot = sns.barplot(x= "desc_x", y = "mean", data =arch.sort_values(by = "desc_x"), hue = "arch",
               hue_order = order, palette =arch_palette,) #showfliers = False )

for n, p in enumerate(splot.patches):
    value = counts.iloc[n]["count"].astype(int)
    splot.annotate(value,
                   (p.get_x() + p.get_width() / 2.,0.01),
                   ha = 'center', va = 'baseline',
                   size=15,
                   rotation = 90,
                   color = "white",
                   xytext = (0, 1),
                   textcoords = 'offset points'
                   )
ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
ax.set_ylabel("Mean LINSIGHT score")
ax.set_xlabel("")
ax.legend(bbox_to_anchor = (1,1))
plt.savefig("%sFigS32_ROADMAP_0.1sampled_mean_linsight_overall_no_exon_relative.pdf"%RE, bbox_inches = "tight")


# In[22]:


arch.groupby('arch')["mean"].mean()


# In[ ]:
