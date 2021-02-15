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



# In[45]:


def plot_reg(taxon2, df, n_bins):

    test = df.loc[df.taxon2 == taxon2] # get age-specfic dataframe

    fig, ( ax2) = plt.subplots(figsize = (8,8))
    sns.set("talk")
    # architecture-specific dataframe
    simple = test.loc[test.arch == "simple"]

    complexenh = test.loc[test.arch != "simple"]

    x = "tissue_overlap"
    y = "linsight_score"


    # do linear regression for simple enhancers in age - pleiotropy x length
    # get coeffs of linear fit

    slope, intercept, r_value, p_value, std_err = stats.linregress(simple[x],simple[y])


    # get counts
    counts = df.loc[df.taxon2 == taxon2,["arch", "enh_id", "taxon2"]].drop_duplicates()
    counts = counts.groupby("arch")["enh_id"].count().reset_index()
    nsimple = counts.loc[counts.arch == "simple", "enh_id"].iloc[0]
    ncomplex = counts.loc[counts.arch != "simple", "enh_id"].iloc[0]


    # plot regplot w/ linear regression annotation
    sns.regplot(x=x, y=y, data = simple, x_bins = n_bins,fit_reg=True,
                line_kws={'label':"r2=%s" % round(r_value, 2)},
                truncate = False,
                #line_kws={'label':"r2={0:.3f}".format(r_value)},
                #line_kws={'label':"y={0:.3f}x+{1:.3f}".format(slope,intercept)},
                color="y", ax = ax2,
                x_estimator=np.mean)


    # plot complex enhancers pleiotropy x length
    if taxon2 != "Homo sapiens (0)":

        slope, intercept, r_value, p_value, std_err = stats.linregress(complexenh[x],complexenh[y])


        sns.regplot(x=x, y=y, data = complexenh, x_bins = n_bins,fit_reg=True,
        truncate = False,
                line_kws={'label':"r2=%s" % round(r_value, 2)},
                #line_kws={'label':"y={0:.3f}x+{1:.3f}".format(slope,intercept)},
                color="g", ax = ax2,
                x_estimator=np.mean)

        ax2.set(title = "%s" % taxon2,
                #xlim = (0,complexenh.tissue_overlap.max()),
                xlabel = 'context overlap\n nsimple=%s, ncomplex=%s' % (nsimple, ncomplex),
                ylabel = "LINSIGHT")

        #ax2.set_xticks(np.arange(0, max.loc[max.tissue_overlap.max(), step = 10))
        ax2.legend()
    plt.tight_layout()
    return fig


# In[46]:


taxon2 = "Tetrapoda (352)"
nbins = 10
plot_reg(taxon2, df, nbins)
#%%

df.sort_values(by="mrca_2")

for taxon2 in df.sort_values(by="mrca_2").taxon2.unique():
    print(taxon2)
    if taxon2 != "Homo sapiens (0)" and taxon2 != -1:
        nbins = 10
        fig = plot_reg(taxon2, df, nbins)
        sid = taxon2.split(" ")[0]
        plt.savefig("%sfantom_trimmedpleiotropy_x_raw_LINSIGHT_%s.pdf" %(RE, sid))
