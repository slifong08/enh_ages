import glob
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.ticker as ticker
import numpy as np
import os, sys
from scipy import stats
import seaborn as sns
import statsmodels
import statsmodels.api as sm

RE ="/dors/capra_lab/projects/enhancer_ages/fantom/results/for_publication/age_breaks/"

es_colors = [ "slate grey","greyish"]
es_pal = sns.xkcd_palette(es_colors)
sns.palplot(es_pal)

#%% Files
path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/"

enh = "%sFANTOM_enh_age_arch_full_matrix.tsv" % path
summaryEnh = "%sFANTOM_enh_age_arch_summary_matrix.tsv" % path

shuf = "%sSHUFFLED_FANTOM_enh_age_arch_full_matrix.tsv" % path
summaryShuf = "%sSHUFFLE_FANTOM_enh_age_arch_summary_matrix.tsv" % path

#%% other summary files

# age and taxon file
syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2", "mya", "mya2"]] # whittle down the df
syn_gen_bkgd

# tissue/cell line dataset descriptions
desc_file = "/dors/capra_lab/data/fantom/fantom5/facet_expressed_enhancers/sample_id_descriptions.txt"
desc_df= pd.read_csv(desc_file, sep = '\t', header = None)

#%% LOAD Files


shuffle = pd.read_csv(shuf, sep = '\t')
final_merge = pd.read_csv(enh, sep = '\t')


#%%


enh_lens = final_merge.groupby(["enh_id", "core_remodeling","enh_len", "arch"])["mrca_2"].max().reset_index()
enh_lens.mrca_2 = enh_lens.mrca_2.round(3)

enh_lens = pd.merge(enh_lens, syn_gen_bkgd[[ "mrca_2","taxon2", "mya2"]], how = "left", on = "mrca_2")
enh_lens = enh_lens.drop_duplicates()
print(enh_lens.shape)

enh_lens["datatype"]="FANTOM"
enh_lens.head()


#%%


shuf_len = shuffle.groupby(["enh_id", "core_remodeling",\
 "enh_len", "shuf_id", "arch"])["mrca_2"].max().reset_index()

print(shuf_len.shape)

shuf_len["datatype"]="SHUFFLE"
shuf_len.mrca_2 = shuf_len.mrca_2.round(3)
shuf_len = pd.merge(shuf_len, syn_gen_bkgd[["mrca_2", "taxon2", "mya2"]],\
 how = "left", on = "mrca_2")

shuf_len = shuf_len.drop_duplicates()
print(shuf_len.shape)


#%% sample 1/10th of shuffle dataset and combine w/ enhancers for comparison


sample_shuf = shuf_len.sample(frac = 0.1)

lens = pd.concat([enh_lens, sample_shuf])


#%% calculate frequency of obs and expected enhancer lengths from fantom, shuffle
x = len(enh_lens) # how many enhancers are there?
enh_len_freq = enh_lens.groupby(["mrca_2", "core_remodeling"])["enh_id"]\
.count().divide(x).round(3).reset_index()
enh_len_freq.columns = ["mrca_2", "core_remodeling", "enh_freq"]

y = len(shuf_len) # how many shuffled enhancers are there?
shuf_len_freq = shuf_len.groupby(["mrca_2", "core_remodeling"])["enh_id"]\
.count().divide(y).round(3).reset_index()
shuf_len_freq.columns = ["mrca_2", "core_remodeling", "shuf_freq"]

# calculat difference between exp and observed frequency
len_freq = pd.merge(enh_len_freq, shuf_len_freq, how = 'left' )
len_freq["dif"] = len_freq.enh_freq / len_freq.shuf_freq
len_freq

shuf_len = shuf_len.loc[~(shuf_len.arch.str.contains("complex")&(shuf_len.mrca_2 == 0.0))]
enh_lens = pd.merge(enh_lens, len_freq, how = "left")
shuf_len = pd.merge(shuf_len, len_freq, how = "left")
shuf_len.head()

#%% plot enhancer architecture length per age
e_colors = [ "amber", "faded green"]
e_pal = sns.xkcd_palette(e_colors)
s_colors = [ "greyish", "slate grey"]
s_pal = sns.xkcd_palette(s_colors)

hue_order = ["FANTOM", "Shuffle"]
fig,(ax1) = plt.subplots(figsize = (8, 8))
order =["Simple", "Complexenh"]
sns.barplot(y = "enh_len", x = "taxon2",\
data = enh_lens.sort_values(by = "mrca_2"), ax = ax1,\
hue = "arch",  palette = e_pal, estimator = np.median)#showfliers=False)


ms, msp = stats.mannwhitneyu(enh_lens.enh_len.loc[enh_lens.arch.str.contains("imple")],
                            shuf_len.enh_len.loc[shuf_len.arch.str.contains("imple")])
print("simple", ms, msp)

mc, mcp = stats.mannwhitneyu(enh_lens.enh_len.loc[enh_lens.arch.str.contains("omplex")],
                            shuf_len.enh_len.loc[shuf_len.arch.str.contains("omplex")])
print("complex", mc, mcp)
ax1.set(ylabel= "Enhancer Length (bp)", ylim = (190,400), xlabel = "")
ax1.set_xticklabels(ax1.get_xticklabels(), rotation = 90, horizontalalignment = "left")


ax1.get_legend().remove()
plt.savefig("%sfig2c-Fantom_ENH_MRCA_x_LEN_ENH.pdf" % RE, bbox_inches = "tight")

""" RESULTS enhancer lengths v. expected shuffle lengths for simple, complex
simple 15249058638.5 2.627822101951775e-05
complex 6146319551.5 1.3788416832439736e-06
"""
#%% Supplemental Figure 2.6 - FANTOM/SHUFFLE enhancer age x arch length
len_concat = pd.concat([enh_lens, shuf_len])

len_concat["datatype2"] = len_concat.datatype + "-"+ len_concat.arch
len_concat["datatype2"].unique()


e_colors = [ "amber","greyish", "faded green", "slate grey"]
e_pal = sns.xkcd_palette(e_colors)
hue_order = ['FANTOM-simple','SHUFFLE-simple','FANTOM-complexenh','SHUFFLE-complexenh']

#%%
fig,(ax1) = plt.subplots(figsize = (8, 8))
order =["Simple", "Complexenh"]
sns.barplot(y = "enh_len", x = "taxon2",
            data = len_concat.sort_values(by = "mrca_2"), ax = ax1,
            hue = "datatype2",  palette = e_pal,
            hue_order = hue_order,
            estimator = np.median)#showfliers=False)


ax1.set(ylabel = "Enhancer Length (bp)", xlabel = "")
ax1.set_xticklabels(ax1.get_xticklabels(), rotation = 90, horizontalalignment = "left")

ax1.legend().remove()
plt.savefig("%sfigS2.6A-Fantom_ENH_MRCA_x_LEN_ENH.pdf" % RE, bbox_inches = "tight")
#%%
e_colors = [ "amber","greyish", "faded green", "slate grey"]
e_pal = sns.xkcd_palette(e_colors)
hue_order = ['FANTOM-simple','SHUFFLE-simple','FANTOM-complexenh','SHUFFLE-complexenh']
#fig,(ax1) = plt.subplots(figsize = (8, 8))
order =["Simple", "Complexenh"]
sns.lmplot(y = "enh_len", x = "mya2", data = len_concat,\
            hue = "datatype2",  palette = e_pal,
            hue_order = hue_order, scatter = False)
plt.savefig("%sfigS2.6B-LM_enh_len.pdf" % RE, bbox_inches = "tight")
#%%
