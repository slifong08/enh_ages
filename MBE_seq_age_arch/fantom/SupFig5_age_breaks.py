import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

RE ="/dors/capra_lab/projects/enhancer_ages/fantom/results/for_publication/age_breaks/"
cs = ["faded green", "greyish"]
cs_pal = sns.xkcd_palette(cs)
sns.palplot(cs_pal)
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

# ENHANCER DATAFRAME
core_breaks = final_merge.groupby(["enh_id", "core_remodeling"])\
["mrca_2", "seg_index"].max().reset_index()

core_breaks = pd.merge(core_breaks, syn_gen_bkgd[["mrca_2", "taxon2"]])

core_breaks["datatype"] = "FANTOM"

core_breaks = core_breaks.drop_duplicates()

# SHUFFLE DATAFRAME
smrca = shuffle.groupby(["enh_id", "core_remodeling", "shuf_id"])\
["mrca_2"].max().reset_index()

# DATA CLEAN UP - remove the enh_id errors in the shuffle db
# where the oldest complexenh age is homosapiens
smrca_list = list(smrca.loc[(smrca.core_remodeling == 1) \
                            & (smrca.mrca_2 == 0.00)]["enh_id"])

smrca = smrca[~smrca.enh_id.isin(smrca_list)]

sseg = shuffle.groupby(["enh_id", "core_remodeling", "shuf_id"])\
["seg_index"].max().reset_index()

shuf_core_breaks = pd.merge(smrca, sseg, how ="left",\
on = ("enh_id", "core_remodeling", "shuf_id"))

shuf_core_breaks = pd.merge(shuf_core_breaks,\
 syn_gen_bkgd[["mrca_2", "taxon2"]])
shuf_core_breaks["datatype"] = "SHUFFLE"
shuf_core_breaks=shuf_core_breaks.drop_duplicates()

plot_breaks = pd.concat([core_breaks,shuf_core_breaks])


#%% plot break means for enhancers, shuffle
f, ax = plt.subplots(figsize = (8,8))
sns.set_context("poster")

hue_order = ["FANTOM", "SHUFFLE"]

sns.barplot(x = "datatype", y = "seg_index",
            data =plot_breaks[plot_breaks.core_remodeling ==1],
            palette = cs_pal)

ax.set_xticklabels(ax.get_xticklabels(), rotation = 90, horizontalalignment = "left")

ax.set(ylabel = "# age segments", title = "FANTOM age segments v MRCA",
 xlabel = "", ylim = (0,3.2))

ax.legend().remove()

plt.savefig("%sFigS2.4A-FANTOM_age_segments_bar.pdf" %RE, bbox_inches = "tight")


#%% plot age x break means for enhancers, shuffle

f, ax = plt.subplots(figsize = (8,8))

hue_order = ["FANTOM", "SHUFFLE"]
sns.barplot(x = "taxon2", y = "seg_index",
            data =plot_breaks[plot_breaks.core_remodeling ==1].sort_values(by = "mrca_2"),
            palette = cs_pal, hue = "datatype", hue_order = hue_order)

ax.set_xticklabels(ax.get_xticklabels(), rotation = 90, horizontalalignment = "left")
ax.set(ylabel = "# age segments", title = "FANTOM age segments v MRCA",
xlabel = "", ylim = (0,3.2))
ax.get_legend().remove()

plt.savefig("%sFigS2.4B-FANTOM_age_segments_v_MRCA_bar.pdf" %RE, bbox_inches = "tight")
#%% PERFORM A BUNCH OF KS tests
from scipy.stats import mstats
kw_list = []
fantom_list = []
shuffle_list = []
for i in plot_breaks.mrca_2.unique():
    for d in plot_breaks.datatype.unique():
        mrca_list = plot_breaks.loc[(plot_breaks.mrca_2 == i)
                                    & (plot_breaks.datatype == d)
                                    & (plot_breaks.core_remodeling ==1) , "seg_index"].to_list() # collect the # segments per age
        kw_list.append(mrca_list)
        if d == "FANTOM":
            fantom_list.append(mrca_list) # collect the # segments per age fantom
        elif d == "SHUFFLE":
            shuffle_list.append(mrca_list)# collect the # segments per age shuffle

args=[l for l in kw_list]
argsF=[l for l in fantom_list]
argsS=[l for l in shuffle_list]

#%% KW
# KW on age segment differences for all age and datatype complex enhancers
print(stats.mstats.kruskalwallis(*args))

# KW on age segment differences among FANTOM complex enhancers
print(stats.mstats.kruskalwallis(*argsF))
# KW on age segment differences among Shuffle complex enhancers
print(stats.mstats.kruskalwallis(*argsS))
"""RESULTS
KW among all age and datatype complex enhancers
KruskalResult(statistic=42210.25967748665, pvalue=0.0)

KW among all FANTOM complex enhancers ages x breaks
KruskalResult(statistic=436.8441208479754, pvalue=1.8594342780202455e-88)

KW among all SHUFFLE complex enhancers ages x breaks
KruskalResult(statistic=41595.42171022936, pvalue=0.0)
"""
