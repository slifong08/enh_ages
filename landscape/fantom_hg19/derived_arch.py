import matplotlib.pyplot as plt
import os
import pandas as pd
from scipy import stats
import seaborn as sns
import statsmodels
import statsmodels.api as sm
import subprocess

colors = [ "amber", "faded green"]
palette = sns.xkcd_palette(colors)
sns.palplot(palette)

colors = [ "dusty blue", "greyish"]
es = sns.xkcd_palette(colors)
sns.palplot(es)

colors = [ "dusty purple", "grey"]
pur = sns.xkcd_palette(colors)
sns.palplot(pur)


FANTOMPATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/all_fantom_enh/ages/"
FANTOMFILE = "syn_breaks_all_fantom_enh_ages.bed"
FANTOM = os.path.join(FANTOMPATH, FANTOMFILE)

SHUFPATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/shuffle/first_round_breaks"
SHUFFILE = "no-exon_syn_breaks_shuf-all_fantom_enh_ages.bed"
SHUFF = os.path.join(SHUFPATH, SHUFFILE)


RE = "/dors/capra_lab/projects/enhancer_ages/landscape/results/fantom/"


#%%


def format_syndf(enh_age_file):

    syn_cols = ["chr_syn", "start_syn", "end_syn",
    "enh_id",
    "chr", "start", "end",
    "seg_index", "core_remodeling", "core",
    "mrca",]

    syn = pd.read_csv(enh_age_file, sep ='\t', header = None, names = syn_cols)

    syn["syn_id"] = syn.chr_syn + ":" + syn.start_syn.map(str) + "-" + syn.end_syn.map(str)

    syn["syn_len"] = syn.end_syn - syn.start_syn
    syn["enh_len"] = syn.end - syn.start


    # age and taxon file
    syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
    syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
    syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

    syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2"]] # whittle down the df
    syn["mrca"] = syn["mrca"].round(3) # round the ages

    syn = pd.merge(syn, syn_gen_bkgd, how = "left", on = "mrca")

    syn["arch"] = "simple"

    syn.loc[(syn.core_remodeling ==1) & (syn.core ==1), "arch"] = "complex_core"
    syn.loc[(syn.core_remodeling ==1) & (syn.core ==0), "arch"] = "complex_derived"

    return syn

#%%


enh = format_syndf(FANTOM)
enh["id"] = "FANTOM"

shuf = format_syndf(SHUFF)
shuf["id"] = "SHUFFLE"
shuf.head()

df = pd.concat([enh, shuf])
df.shape


#%% make df of only cores


core = df.loc[df.core ==1].copy()
core.head()


#%%


core["pct_core"] = core.syn_len.divide(core.enh_len)
core["pct_der"] = 1 - core["pct_core"]

#%%


x = "arch"
y = "pct_der"
hue = "id"
data = core.loc[core.arch == "complex_core"]

sns.set("poster")

fig, ax = plt.subplots(figsize=(6,6))

# plot
sns.barplot(x = x, y = y, data = data, hue = hue, palette =es, errwidth= 10)

ax.set(ylabel = "Percent of enhancer length", xlabel = "")
ax.set_xticklabels(["derived"])
ax.legend(bbox_to_anchor = (1,1))
plt.savefig("%sfantom_percent_der.pdf" % RE, bbox_inches = "tight")

#%%


x = "arch"
y = "pct_core"
hue = "id"
data = core

sns.set("poster")

fig, ax = plt.subplots(figsize=(6,6))
# plot cores

sns.barplot(x = x, y = y, data = data, hue = hue, palette =pur, errwidth= 10)

ax.set(ylabel = "Percent of enhancer length", xlabel = "")
ax.legend(bbox_to_anchor = (1,1))
plt.savefig("%sfantom_percent_core_only.pdf" % RE, bbox_inches = "tight")

#%%
"""
#core.groupby(["arch", "id"])["pct_core"].median()

arch          id
complex_core  FANTOM     0.485368
              SHUFFLE    0.398585
simple        FANTOM     1.000000
              SHUFFLE    1.000000


"""

#%%
x = "mrca_2"
y = "pct_core"
hue = "id"
data = core.loc[core.arch == "complex_core"]

sns.set("poster")

xlabs = ["prim", "euar", "bore", "euth", "ther", "mam", "amni", "tetr", "vert"]
fig, ax = plt.subplots(figsize=(6,6))
sns.barplot(x = x, y = y, data = data, hue = hue, palette = pur, errwidth= 5)
ax.set(ylabel = "Percent of enhancer length", xlabel = "")

ax.set_xticklabels(xlabs, rotation = 90)
ax.legend(bbox_to_anchor = (1,1))
plt.savefig("%sfantom_percent_core_mrca_2.pdf" % RE, bbox_inches = "tight")

#%%

x = "mrca_2"
y = "pct_der"
hue = "id"
data = core.loc[core.arch == "complex_core"]

sns.set("poster")

xlabs = ["prim", "euar", "bore", "euth", "ther", "mam", "amni", "tetr", "vert"]
fig, ax = plt.subplots(figsize=(6,6))
sns.barplot(x = x, y = y, data = data, hue = hue, palette = es, errwidth= 5)
ax.set(ylabel = "Percent of enhancer length", xlabel = "")
ax.set_xticklabels(xlabs, rotation = 90)
ax.legend(bbox_to_anchor = (1,1))
plt.savefig("%sfantom_percent_der_mrca_2.pdf" % RE, bbox_inches = "tight")
