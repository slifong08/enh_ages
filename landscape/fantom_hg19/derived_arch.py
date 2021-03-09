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

colors = [ "amber", "greyish", "faded green", "grey"]
enhpal = sns.xkcd_palette(colors)
sns.palplot(enhpal)

colors = [ "amber", "greyish",  "dusty purple", "brown grey",  "windows blue", "bluey grey"]
archpal = sns.xkcd_palette(colors)
sns.palplot(archpal)

FANTOMPATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/all_fantom_enh/ages/"
FANTOMFILE = "syn_breaks_all_fantom_enh_ages.bed"
FANTOM = os.path.join(FANTOMPATH, FANTOMFILE)

SHUFPATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/shuffle/first_round_breaks"
SHUFFILE = "no-exon_syn_breaks_shuf-all_fantom_enh_ages.bed"
SHUFF = os.path.join(SHUFPATH, SHUFFILE)


RE = "/dors/capra_lab/projects/enhancer_ages/landscape/results/fantom/"


#%%


def add_arch_labels(df):

    df["arch"] = "complex_core"

    df.loc[df.core_remodeling ==0, "arch"] = "simple"
    df.loc[df.core ==0, "arch"] = "complex_derived"

    return df


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

    labeled_syn = add_arch_labels(syn) # add architecture labels

    return labeled_syn


def clean_shuffles(df):

    remove_list = []
    shuf_remove_ids = df.loc[(df["pct_enh"] ==1) & (df.core ==0)]

    df = df.loc[~df.enh_id.isin(shuf_remove_ids)]

    return df


#%%


enh = format_syndf(FANTOM)
enh["id"] = "FANTOM"


shuf = format_syndf(SHUFF)

shuf["id"] = "SHUFFLE"
shuf.head()

df = pd.concat([enh, shuf])
df.shape



#%% summarize architecture lengths per enhancer


sum_archlen = df.groupby(['id', 'enh_id', 'enh_len', 'core_remodeling', 'core'])["syn_len"].sum().reset_index().drop_duplicates()

sum_archlen = add_arch_labels(sum_archlen) # add architecture labels

 # calculate percent arch per enhancer
sum_archlen["pct_enh"] = sum_archlen.syn_len.divide(sum_archlen.enh_len)
sum_archlen = clean_shuffles(sum_archlen, shuf) # remove the shuffles w/o core

# add oldest ages of enhancer information
enh_mrcas = df.groupby("enh_id")[["mrca_2", "taxon2"]].max().reset_index()
sum_archlen = pd.merge(sum_archlen, enh_mrcas, how = "left", on = "enh_id")

sum_archlen.head()

sum_archlen = sum_archlen[~sum_archlen["enh_id"].isin(shuf_remove_ids)]
#%%

sum_archlen["arch_"] = "simple"
sum_archlen.loc[sum_archlen.core_remodeling ==1, "arch_"] = "complex"
sum_archlen["arch__"] = sum_archlen["arch_"] + "-" + sum_archlen['id']
sum_archlen["arch__"].unique()

#%% plot syn lengths per age, arch.


x = "mrca_2"
y = "syn_len"
data = sum_archlen.sample(frac = 0.25)
hue = "arch__"


xlabs = ["homo", "prim", "euar", "bore", "euth", "ther", "mam", "amni", "tetr", "vert"]
order = ["simple-FANTOM", 'simple-SHUFFLE', 'complex-FANTOM', 'complex-SHUFFLE']

fig, ax = plt.subplots(figsize = (12,6))
sns.barplot(x=x, y=y, data = data,
            hue = hue,
            palette = enhpal, hue_order = order )

ax.set(ylabel = "syntenic length", xlabel = "")
ax.set_xticklabels(xlabs, rotation = 90)
ax.legend(bbox_to_anchor = (1,1))
outf = "%smrca_x_syn_lengths_arch.pdf" % RE

plt.savefig(outf, bbox_inches = 'tight')


#%% plot length of segments


ticklabs = ["simple", "complex\ncore", "complex\nderived"]

x = "arch"
y = "syn_len"
hue = "id"
data = sum_archlen
order = ["simple", "complex_core", "complex_derived"]

fig, ax = plt.subplots(figsize = (6,6))

sns.barplot(x=x, y=y, data = data, hue = hue, palette = es, order = order)
ax.set_xticklabels(ticklabs)
ax.legend(bbox_to_anchor = (1,1))
ax.set(ylabel = "sum length (bp)", xlabel = "")

plt.savefig("%ssum_arch_length_all_fantom_shuffle.pdf" %RE, bbox_inches = "tight")


#%% mean lengths


sum_archlen.groupby(["id", "arch"])["syn_len"].mean()
'''
mean lengths

id       arch
FANTOM   complex_core       17.824845 bp long
         complex_derived    173.842735
         simple             276.605902

SHUFFLE  complex_core       177.303235
         complex_derived    189.748677
         simple             274.325146
'''

sum_archlen.groupby(["id", "arch"])["pct_enh"].mean()
#%%



x = "arch"
y = "pct_enh"
hue = "id"


sns.set("poster")

fig, ax = plt.subplots(figsize=(6,6))

# plot
sns.barplot(x = x, y = y, data = data, order = order,
hue = hue, palette =archpal, errwidth= 10)

ax.set(ylabel = "Percent of enhancer length", xlabel = "")
ax.set_xticklabels(ticklabs)
ax.legend(bbox_to_anchor = (1,1))
plt.savefig("%sfantom_percent_all.pdf" % RE, bbox_inches = "tight")




#%%
sum_archlen["id2"] = sum_archlen["arch"] + "-" + sum_archlen["id"]


x = "mrca_2"
y = "pct_enh"
hue = "id2"
order2 = ["simple-FANTOM", "simple-SHUFFLE",
"complex_core-FANTOM","complex_core-SHUFFLE",
 "complex_derived-FANTOM", "complex_derived-SHUFFLE"]

sns.set("poster")


fig, ax = plt.subplots(figsize=(12,6))
sns.barplot(x = x, y = y, data = data, hue = hue, hue_order = order2,
palette = archpal, errwidth= 5)

ax.set(ylabel = "Percent of enhancer length", xlabel = "")

ax.set_xticklabels(xlabs, rotation = 90)
ax.legend(bbox_to_anchor = (1,1))
plt.savefig("%sfantom_percent_arch_mrca_2.pdf" % RE, bbox_inches = "tight")

#%%

x = "mrca_2"
y = "syn_len"
hue = "id2"


sns.set("poster")


fig, ax = plt.subplots(figsize=(12,6))
sns.barplot(x = x, y = y, data = data, hue = hue,
hue_order = order2,  palette = archpal, errwidth= 5)

ax.set(ylabel = "syntenic length", xlabel = "")
ax.set_xticklabels(xlabs, rotation = 90)
ax.legend(bbox_to_anchor = (1,1))
plt.savefig("%sfantom_syn_len_mrca_2.pdf" % RE, bbox_inches = "tight")
#%%

shuf_remove_ids = sum_archlen.loc[(sum_archlen["mrca_2"] == 0) & (sum_archlen.core ==0), "enh_id"]

test = "chr14:89971043-89971729"
sum_archlen.loc[sum_archlen.enh_id == test]
shuf.loc[shuf.enh_id == test]
testsyn = "chr14:89971596-89971729"
shuf.loc[shuf.syn_id == testsyn]
