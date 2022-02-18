
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

RE ="/dors/capra_lab/projects/enhancer_ages/fantom/results/for_publication/syn/"

arch_colors = ["amber", "dusty purple", "windows blue","greyish"]
arch_palette = sns.xkcd_palette(arch_colors)
sns.palplot(arch_palette)

cs = ["faded green", "greyish"]
cs_pal = sns.xkcd_palette(cs)
sns.palplot(cs_pal)

fullpal_colors = [ "dusty purple", "slate grey", "windows blue", "blue grey", "amber", "greyish"]
full_palette = sns.xkcd_palette(fullpal_colors)
sns.palplot(full_palette)

fullenh_colors = [ "faded green", "slate grey", "amber", "greyish"]
fullenh_palette = sns.xkcd_palette(fullenh_colors)
sns.palplot(fullenh_palette)

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
enh = "%sFANTOM_enh_age_arch_full_matrix.tsv" % path
summaryEnh = "%sFANTOM_enh_age_arch_summary_matrix.tsv" % path

shuf = "%sSHUFFLED_FANTOM_enh_age_arch_full_matrix.tsv" % path
summaryShuf = "%sSHUFFLE_FANTOM_enh_age_arch_summary_matrix.tsv" % path

shuffle = pd.read_csv(shuf, sep = '\t')
shuffle.mrca_2 = shuffle.mrca_2.round(3)

final_merge = pd.read_csv(enh, sep = '\t')
final_merge.mrca_2 = final_merge.mrca_2.round(3)


#%%
final_merge.head()
#%% enhancer + shuffle syntenic length analysis

# Shuffled syntenic lengths
ssyn_lens = shuffle[["syn_id", "code","syn_len", "mrca_2"]].drop_duplicates()
ssyn_lens["dataset"] = "shuffle"

syn_lens = final_merge[["syn_id", "code","syn_len", "mrca_2"]].drop_duplicates()
syn_lens["dataset"] = "fantom"
plot = pd.concat([syn_lens, ssyn_lens])


#%%

order = ["simple", "complex_core", "derived"]
hue_order = ["fantom", "shuffle"]
fig, ax = plt.subplots(figsize = (8,8))
sns.set_context("poster")
sns.boxplot(x = "code", y = "syn_len",
            data = plot, palette = cs_pal, order = order,
            hue = "dataset",
            hue_order = hue_order,
           showfliers = False)

ax.set(xlabel = "Architecture", ylabel= "Enhancer Length (bp)",
xticklabels = ["Simple", "Complex\nCore", "Derived"])

plt.savefig("%sfantom_v_shuffle_syn_lens.pdf" %RE, bbox_inches = "tight")

syn_lens_df = syn_lens.groupby(["code"])["syn_len"].median().reset_index()
print("fantom", syn_lens_df)

ssyn_lens_df = ssyn_lens.groupby(["code"])["syn_len"].median().reset_index()
print("shuffle", ssyn_lens_df)
""" RESULTS median lengths
fantom     code  syn_len
0  complex_core      157
1       derived       85
2        simple      259

SHUFFLE    code  syn_len
0  complex_core      129
1       derived       90
2        simple      255
"""

#%% enhancer syntenic lengths only

order = ["simple", "complex_core", "derived"]
fig, ax = plt.subplots(figsize = (8,8))

sns.boxplot(x = "code", y = "syn_len",
            data = syn_lens, palette = arch_palette, order = order,
           showfliers = False)
#ax.set_xticklabels(["Simple", "Complexenh"])
ax.set(xlabel = "Architecture", ylabel= "Enhancer Length (bp)",
xticklabels = ["Simple", "Complex\nCore", "Derived"])

ax.set_xticklabels(["Simple", "Complex\nCore", "Derived"])
plt.savefig("%sfantom_syn_lens.pdf" %RE, bbox_inches = "tight")



#%%
plot["dataset2"] = plot.code + "-" + plot.dataset
order = ["simple", "complex_core", "derived"]

fig, ax = plt.subplots(figsize = (16,8))
sns.set_context("poster")
sns.barplot(x = "mrca_2", y = "syn_len",
            data = plot.sort_values(by = "dataset2"),
            hue = "dataset2",
           #showfliers = False,
           palette = full_palette)

ax.set(xlabel = "Architecture", ylabel= "Enhancer Length (bp)")
ax.legend(bbox_to_anchor = (1,1))

plt.savefig("%sfantom_v_shuffle_syn_lens_per_mrca.pdf" %RE, bbox_inches = "tight")

syn_lens_df = syn_lens.groupby(["code"])["syn_len"].median().reset_index()
print("fantom", syn_lens_df)

ssyn_lens_df = ssyn_lens.groupby(["code"])["syn_len"].median().reset_index()
print("shuffle", ssyn_lens_df)


#%%


ssyn_lens_arch = shuffle[["syn_id", "arch","syn_len", "mrca_2"]].drop_duplicates()
ssyn_lens_arch["dataset"] = "shuffle"

syn_lens_arch = final_merge[["syn_id", "arch","syn_len", "mrca_2"]].drop_duplicates()
syn_lens_arch["dataset"] = "fantom"
plot_arch = pd.concat([syn_lens_arch, ssyn_lens_arch])


#%%

plot_arch.arch.unique()

#%%

order = ["simple", "complexenh",]
hue_order = ["fantom", "shuffle"]
fig, ax = plt.subplots(figsize = (8,8))
sns.set_context("poster")
sns.barplot(x = "arch", y = "syn_len",
            data = plot_arch, palette = cs_pal, order = order,
            hue = "dataset",
            #showfliers = False,
            hue_order = hue_order,)

ax.set(xlabel = "Architecture", ylabel= "Enhancer Length (bp)",
xticklabels = ["Simple", "Complex\nenh"])
plt.savefig("%sfantom_v_shuffle_enh_lens_arch.pdf" %RE, bbox_inches = "tight")

syn_lens_df_ = syn_lens_arch.groupby(["arch"])["syn_len"].median().reset_index()
print("fantom", syn_lens_df_)

ssyn_lens_df_ = ssyn_lens_arch.groupby(["arch"])["syn_len"].median().reset_index()
print("shuffle", ssyn_lens_df_)

# MWU analysis of simple enh v. shuffle lengths
simpleM, simpleP = stats.mannwhitneyu(syn_lens_arch.loc[(syn_lens_arch.arch == "simple") & \
(syn_lens_arch.dataset == "fantom"), "syn_len"], \
syn_lens_arch.loc[(syn_lens_arch.arch == "simple") & \
(syn_lens_arch.dataset != "fantom"), "syn_len"])

# MWU analysis of complex enh v. shuffle lengths
complexM, complexP = stats.mannwhitneyu(syn_lens_arch.loc[(syn_lens_arch.arch == "complexenh") & \
(syn_lens_arch.dataset == "fantom"), "syn_len"], \
syn_lens_arch.loc[(syn_lens_arch.arch == "complexenh") & \
(syn_lens_arch.dataset != "fantom"), "syn_len"])

print(simpleM, simpleP, complexM, complexP)


#%%


plot_arch["dataset2"] = plot_arch.arch + "-" + plot_arch.dataset
plot_arch = pd.merge(plot_arch, syn_gen_bkgd[["mrca_2", "taxon2"]])
#%%
fig, ax = plt.subplots(figsize = (16,8))
sns.set_context("poster")
sns.barplot(x = "taxon2", y = "syn_len",
            data = plot_arch.sort_values(by = ["mrca_2", "dataset2"]),
            hue = "dataset2",
           #showfliers = False,
           palette = fullenh_palette)

ax.set(xlabel = "Architecture", ylabel= "Enhancer Length (bp)",
 title= "syntenic blocks")
ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
ax.legend(bbox_to_anchor = (1,1))

plt.savefig("%sfantom_v_shuffle_enh_lens_per_mrca.pdf" %RE, bbox_inches = "tight")
