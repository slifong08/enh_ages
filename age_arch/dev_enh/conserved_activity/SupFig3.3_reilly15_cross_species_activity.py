import glob
import matplotlib.pyplot as plt
import os, sys
import pandas as pd
import numpy as np
import seaborn as sns
from scipy import stats

colors = [ "amber", "faded green"]
palette = sns.xkcd_palette(colors)
sns.palplot(palette)
plt.rcParams.update({'font.size': 15})
sns.set_style("white")

#%%


RE = "/dors/capra_lab/projects/enhancer_ages/reilly15/results/pleiotropy/"
path = "/dors/capra_lab/projects/enhancer_ages/reilly15/data/multiintersect/"
f = "%strim-0.5_multiintersect_hu_count.bed" % path


#%%

def custom_round(x, base=10):
    return int(base * round(float(x)/base))


def match_len(simple_df, complex_df, base_len):

    columns = ["enh_id", "enh_len"]
    columns_names = ["matching_ids", "matching_len"]

    simple = simple_df[columns].drop_duplicates()

    simple.columns = columns_names
    simple.matching_len = simple.matching_len.astype(float).apply(lambda x: custom_round(x, base=base_len)) # round to the nearest 100bp

    complex = complex_df[columns].drop_duplicates()
    complex.columns = columns_names

    complex.matching_len = complex.matching_len.astype(float).apply(lambda x: custom_round(x, base=base_len))

    lens = set(list(simple.matching_len.unique()) + list(complex.matching_len.unique()))

    match_dict = {}

    for length in lens:
        complexn = len(complex.loc[complex.matching_len == length])
        simplen = len(simple.loc[simple.matching_len == length])

        sn = min(complexn, simplen)

        if length > 0 and sn > 0:


            # find 100 matched enhancer samples

            complex_ids = complex.loc[complex.matching_len == length].sample(n = complexn, replace = False) # sample w/ replacement
            simple_ids = simple.loc[simple.matching_len == length].sample(n = simplen, replace = False) # sample w/ replacement
            balanced = pd.concat([simple_ids, complex_ids])
            match_dict[sn] = balanced

    final_matched_id = pd.concat(match_dict.values())

    return final_matched_id.matching_ids.unique()

#%%


df = pd.read_csv(f, sep= '\t', header = None)
df.columns = ["chr_enh","start_enh", "end_enh",	"enh_id", "fourth_col", "id",\
"core_remodeling","arch","seg_index","mrca","enh_len", "taxon", "mrca_2",\
"taxon2",	"mya",	"mya2",	"seg_den",	"datatype", "count_overlap"]
df["enh_len"] = df.end_enh - df.start_enh
df["mrca"] = df["mrca"].astype(float).round(3)
df = df.loc[df.chr_enh != "chrX"]
df = df.loc[df.enh_len < 10000]
#df = pd.merge(df, syn_gen_bkgd, how = "left", on = "mrca")
median = df.seg_index.median()
#median = 3
df.loc[df.seg_index >=median, "core_remodeling"] = 1
df.loc[df.seg_index >=median, "arch"] = "complexenh"

df.loc[df.seg_index <median, "core_remodeling"] = 0
df.loc[df.seg_index <median, "arch"] = "simple"

#%% remove human, primate specific sequences. Sequence must be as olds as euarchontaglires (common ancestor with mouse) in order to be evaluated here.
# removes 470 sequences

print(len(df))
print(median)
#%%
df.groupby(["core_remodeling", "arch"])["enh_id"].count()

#%%

simple_df = df.loc[df.core_remodeling ==0]

complex_df = df.loc[(df.core_remodeling ==1) & (df.enh_len<=simple_df.enh_len.max())]

#%%

matched_id = match_len(simple_df, complex_df, 100)


matched = df.loc[df.enh_id.isin(matched_id)]

matched["arch"] = "simple"
matched.loc[matched.core_remodeling ==1, "arch"] = "complexenh"
matched.mrca = matched.mrca.round(3)
#matched = pd.merge(matched, syn_gen_bkgd[["mrca", "mrca_2", "taxon2"]], how = "left", on ="mrca")
matched = matched.loc[matched.chr_enh != "chrX"]
matched.head()
#%%

sns.distplot(matched.loc[matched.core_remodeling == 0, "enh_len"], kde = False, norm_hist = True, label = "simple")
sns.distplot(matched.loc[matched.core_remodeling != 1, "enh_len"], kde = False, norm_hist = True, label = "complex")
plt.legend()



#%% plot simple v. complex
order = ["simple", "complexenh"]
sns.boxplot(x= "arch", y = "count_overlap", data = matched, order = order)

#%%

from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator


fig = plt.figure(figsize = (12, 8))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 3])
ax0 = plt.subplot(gs[0])

sns.barplot(x = "arch", y = "count_overlap", data = matched,
            palette = palette, order = order,
            ax = ax0)
nsimple = len(matched.loc[matched.arch == "simple"])
ncomplex = len(matched.loc[matched.arch != "simple"])
labels = ["n = %s"%nsimple, "n = %s"%ncomplex, ]
ax0.set_xticklabels(labels, rotation = 90)
ax0.set(xlabel="", ylabel ="Number of Active Species", ylim=(0,2))
ax0.yaxis.set_major_locator(MultipleLocator(0.5))

sns.set("poster")

ax2 = plt.subplot(gs[1])
sns.barplot(x = "taxon2", y = "count_overlap", hue = "core_remodeling",
              data = matched.sort_values(by = "mrca_2"),
                palette = palette,
            ax = ax2)

ax2.set_xticklabels(ax2.get_xticklabels(), rotation = 90)
sns.set("poster")
ax2.yaxis.set_major_locator(MultipleLocator(0.5))
ax2.set(ylabel="",  ylim=(0,2))
ax2.legend().remove()
plt.savefig("%sFigS3.3-JOINT_barplot_reilly15_cross_species_overlap_x_mrca_2.pdf" % RE, bbox_inches = "tight" )

#%%
mwu, p = stats.mannwhitneyu(matched.loc[matched.arch == "simple", "count_overlap"],\
matched.loc[matched.arch != "simple", "count_overlap"])
print(mwu, p)


#%%
matched.groupby("arch")["count_overlap"].mean()

#%%
RE
