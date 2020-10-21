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
arch_colors = [ "amber", "dusty purple", "windows blue","greyish"]
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


# In[2]:


linsight_path = "/dors/capra_lab/projects/enhancer_ages/linsight/data/"
fantom_fs = glob.glob("%sFANTOM_*.bed" % linsight_path)

# raw length enh only
fantom_fs = "/dors/capra_lab/projects/enhancer_ages/linsight/data/all_unique_fantom_erna_112_tissue_linsight.bed"
# trimmed length enh only
fantom_fs = "/dors/capra_lab/projects/enhancer_ages/linsight/data/trimmed_all_fantom_enh_112_tissues_multiintersect_0_linsight.bed"


# In[3]:


syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t')
syn_gen_bkgd["mrca"]=syn_gen_bkgd["mrca"].round(3)
syn_gen_bkgd["mrca_2"]=syn_gen_bkgd["mrca_2"].round(3)
syn_gen_bkgd.head()
syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2"]]


# In[39]:



arch_id = (fantom_fs.split("/")[-1]).split(".")[0]
print(arch_id)
df = pd.read_csv(fantom_fs, sep = '\t', header = None, low_memory=False)

df.columns = [ 'chr_enh', 'start_enh','end_enh', "old_len", 'core_remodeling',
              'mrca',"id", "tissue_overlap", "chr_lin",  "start_lin", "end_lin", "linsight_score", "overlap"]


df['enh_len'] = df.end_enh - df.start_enh

df["enh_id"] = df.chr_enh + ":" + df.start_enh.map(str) +"-"+ df.end_enh.map(str)


# arch code

df["arch"] = ""
df.loc[(df.core_remodeling == 0), "arch" ] = "simple"
df.loc[(df.core_remodeling == 1), "arch"] = "complexenh"

df = df.loc[df.arch != ""]

df = df.loc[df.linsight_score != "."] # exclude the loci that do not overlap a linsight score

df.linsight_score = df.linsight_score.astype(float) # turn linsight scores into floats

df["linsight_id"] = df.chr_lin + ":" + df.start_lin.map(str) +"-"+ df.end_lin.map(str)
df.mrca = df.mrca.round(3)
df = pd.merge(df, syn_gen_bkgd, how = "left")

# make a dataframe only based on linsight scores and enhancer architectures

print(len(df.enh_id.unique()))
df.head()


# In[40]:


df.taxon2.unique()


# In[41]:


# take the weighted average for each enhancer
# linsight score * (# overlapping linsight bases/total enhancer length bases)

df['wa'] = df.overlap.divide(df.enh_len)*df.linsight_score

wa = df.groupby(["enh_id", "arch"])["wa"].sum().reset_index() # sum the weighted average per enhancer
wa.wa = wa.wa.round(1)

wa.head()


# In[42]:


sns.boxplot(x ="arch", y = "wa", data = wa, showfliers = False)


# In[43]:


enh_features = df.groupby(["enh_id", "arch"])["mrca_2", "taxon2", "tissue_overlap", "enh_len", "old_len"].max().reset_index()
enh_features.head()


# In[44]:


merged = pd.merge(wa, enh_features, how = "left")
merged.head()


# In[45]:


def plot_reg(taxon2, df):

    test = df.loc[df.taxon2 == taxon2] # get age-specfic dataframe

    fig, ( ax2) = plt.subplots(figsize = (8,8))
    sns.set("talk")
    # architecture-specific dataframe
    simple = test.loc[test.arch == "simple"]

    complexenh = test.loc[test.arch != "simple"]

    x = "tissue_overlap"
    y = "wa"


    # do linear regression for simple enhancers in age - pleiotropy x length
    # get coeffs of linear fit

    slope, intercept, r_value, p_value, std_err = stats.linregress(simple[x],simple[y])


    # plot simple enhancers pleiotropy x length


    # plot regplot w/ linear regression annotation
    sns.regplot(x=x, y=y, data = simple,
                line_kws={'label':"y={0:.3f}x+{1:.3f}".format(slope,intercept)},
                color="y", ax = ax2,
                x_estimator=np.median)


    # plot complex enhancers pleiotropy x length
    if taxon2 != "Homo sapiens (0)":

        slope, intercept, r_value, p_value, std_err = stats.linregress(complexenh[x],complexenh[y])


        sns.regplot(x=x, y=y, data = complexenh,
                line_kws={'label':"y={0:.3f}x+{1:.3f}".format(slope,intercept)},
                color="g", ax = ax2,
                x_estimator=np.median)

        ax2.set(title = "%s" % taxon2,
                xlim = (0,complexenh.tissue_overlap.max()),
                xlabel = 'context overlap',
                ylabel = "LINSIGHT weighted average")

        ax2.set_xticks(np.arange(0, complexenh.tissue_overlap.max(), step = 10))
        ax2.legend()
    plt.tight_layout()
    return fig


# In[46]:


taxon2 = "Tetrapoda (352)"
plot_reg(taxon2, merged)


# In[47]:


merged.taxon2.unique()
merged = merged.sort_values(by="taxon2").fillna(-1)


# In[48]:


merged.sort_values(by="mrca_2")

for taxon2 in merged.sort_values(by="mrca_2").taxon2.unique():
    print(taxon2)
    if taxon2 != "Homo sapiens (0)" and taxon2 != -1:

        fig = plot_reg(taxon2, merged)
        sid = taxon2.split(" ")[0]
        plt.savefig("%sfantom_trimmedpleiotropy_x_LINSIGHT_%s.pdf" %(RE, sid))


#%%
def plot_reg_lin(taxon2, df):

    test = df.loc[df.taxon2 == taxon2] # get age-specfic dataframe

    fig, ( ax2) = plt.subplots(figsize = (8,8))
    sns.set("talk")
    # architecture-specific dataframe
    simple = test.loc[test.arch == "simple"]

    complexenh = test.loc[test.arch != "simple"]

    y = "tissue_overlap"
    x = "wa"


    # do linear regression for simple enhancers in age - pleiotropy x length
    # get coeffs of linear fit

    slope, intercept, r_value, p_value, std_err = stats.linregress(simple[x],simple[y])


    # plot simple enhancers pleiotropy x length


    # plot regplot w/ linear regression annotation
    sns.regplot(x=x, y=y, data = simple,
                line_kws={'label':"y={0:.3f}x+{1:.3f}".format(slope,intercept)},
                color="y", ax = ax2,
                x_estimator=np.median)


    # plot complex enhancers pleiotropy x length
    if taxon2 != "Homo sapiens (0)":

        slope, intercept, r_value, p_value, std_err = stats.linregress(complexenh[x],complexenh[y])


        sns.regplot(x=x, y=y, data = complexenh,
                line_kws={'label':"y={0:.3f}x+{1:.3f}".format(slope,intercept)},
                color="g", ax = ax2,
                x_estimator=np.median)

        ax2.set(title = "%s" % taxon2,
                #xlim = (0,complexenh.tissue_overlap.max()),
                ylabel = 'context overlap',
                xlabel = "LINSIGHT weighted average")

        #ax2.set_xticks(np.arange(0, complexenh.wa.max(), step = 10))
        ax2.legend()
    plt.tight_layout()
    return fig


# In[46]:


taxon2 = "Tetrapoda (352)"
plot_reg_lin(taxon2, merged)


# In[47]:


merged.taxon2.unique()
merged = merged.sort_values(by="taxon2").fillna(-1)


# In[48]:


merged.sort_values(by="mrca_2")

for taxon2 in merged.sort_values(by="mrca_2").taxon2.unique():
    print(taxon2)
    if taxon2 != "Homo sapiens (0)" and taxon2 != -1:

        fig = plot_reg_lin(taxon2, merged)
        sid = taxon2.split(" ")[0]
        plt.savefig("%sfantom_trimmed_LINSIGHT_x_pleiotropy_0.1_%s.pdf" %(RE, sid))