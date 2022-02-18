import glob
import matplotlib.pyplot as plt
import os, sys
import pandas as pd
import seaborn as sns
from scipy import stats

colors = [ "amber", "faded green"]
palette = sns.xkcd_palette(colors)
sns.palplot(palette)
plt.rcParams.update({'font.size': 15})
sns.set_style("white")

#%%


RE = "/dors/capra_lab/projects/enhancer_ages/emera16/results/pleiotropy/"
path = "/dors/capra_lab/projects/enhancer_ages/emera16/data/multiintersect/"
f = "%strim-0.5_multiintersect_hu_count.bed" % path


# age and taxon file
syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2", "mya", "mya2"]] # whittle down the df

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

#chr10	100085975	100088150	chr10:100085975-100088150	chr10:100085975-100088150	3	1	[0, 1, 0]	0.38	1
#%%


df = pd.read_csv(f, sep= '\t', header = None, usecols = [0,1,2,3,5,6,8,9])
df.columns = ["chr_enh","start_enh", "end_enh",	"enh_id",
"seg_index","core_remodeling", "mrca","count_overlap"]

df = df.loc[df["seg_index"] != "max_seg"]
df.seg_index = df.seg_index.astype(int)


#%%

df["enh_len"] = df.end_enh - df.start_enh
df["mrca"] = df["mrca"].astype(float).round(3)
df = df.loc[df.chr_enh != "chrX"]
df = df.loc[df.enh_len < 10000]
df = pd.merge(df, syn_gen_bkgd, how = "left", on = "mrca")

relative_simple = df.seg_index.astype(int).median()
print(relative_simple)

#median = 6
df.loc[df.seg_index >=relative_simple, "core_remodeling"] = 1
df.loc[df.seg_index >=relative_simple, "arch"] = "complexenh"

df.loc[df.seg_index < relative_simple, "core_remodeling"] = 0
df.loc[df.seg_index < relative_simple, "arch"] = "simple"

#%% remove human, primate specific sequences. Sequence must be as olds as euarchontaglires (common ancestor with mouse) in order to be evaluated here.
# removes 470 sequences

print(len(df))

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


xlabs = ["Homo", "Prim", "Euar", "Bore", "Euth", "Ther",
"Mam", "Amni", "Tetr", "Vert"]

# add a blank dataframe so these values get plotted
blankdf = pd.DataFrame({"arch":["simple", "complexenh","simple", "complexenh", "complexenh"],
"count_overlap":[0,0,0,0,0], "core_remodeling":[0,1,0,1,1], "mrca_2":[0, 0, 0.126,  0.126, 0.131],
"taxon2" :[ "Homo sapiens (0)", "Homo sapiens (0)", "Primate (72)", "Primate (72)", "Euarchontoglires (90)"]})

from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator

plot = matched[["arch", "count_overlap", "core_remodeling", "mrca_2", "taxon2", "enh_id"]]

plot = pd.concat([plot, blankdf])

#%%

from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator


fig = plt.figure(figsize = (12, 8))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 3])
ax0 = plt.subplot(gs[0])

splot = sns.barplot(x = "arch", y = "count_overlap", data = plot,
            palette = palette, order = order,
            ax = ax0)
STRAT = 0
agecounts = get_counts(plot, STRAT)

for n, p in enumerate(splot.patches):

    value = agecounts.iloc[n]["enh_id"].astype(int)

    splot.annotate(value,
                   (p.get_x() + p.get_width() / 2.,0.05),
                   ha = 'center', va = 'baseline',
                   size=15,
                   rotation = 90,
                   color = "white",
                   xytext = (0, 1),
                   textcoords = 'offset points'
                   )
ax0.set_xticklabels(["Simple", "Complex"], rotation = 90)
ax0.set(xlabel="", ylabel ="Number of Active Species", ylim=(0,3))
ax0.yaxis.set_major_locator(MultipleLocator(0.5))

sns.set("poster")

ax2 = plt.subplot(gs[1])
mplot = sns.barplot(x = "taxon2", y = "count_overlap", hue = "core_remodeling",
              data = plot.sort_values(by = "mrca_2"),
                palette = palette,
            ax = ax2)

ax2.legend().remove()
STRAT = 1
agecounts = get_counts(plot, STRAT)
agecounts = agecounts.drop_duplicates().reset_index()
agecounts
for n, p in enumerate(mplot.patches):

    value = agecounts.iloc[n]["enh_id"].astype(int)

    mplot.annotate(value,
                   (p.get_x() + p.get_width() / 2.,0.05),
                   ha = 'center', va = 'baseline',
                   size=15,
                   rotation = 90,
                   color = "white",
                   xytext = (0, 1),
                   textcoords = 'offset points'
                   )

ax2.set_xticklabels(ax2.get_xticklabels(), rotation = 90)
sns.set("poster")
ax2.yaxis.set_major_locator(MultipleLocator(0.5))
ax2.set(ylabel="",  ylim=(0,3))

plt.savefig("%sFigS14-JOINT_barplot_emera16_cross_species_overlap_x_mrca_2.pdf" % RE, bbox_inches = "tight" )
#%%
mwu, p = stats.mannwhitneyu(matched.loc[matched.arch == "simple", "count_overlap"],\
matched.loc[matched.arch != "simple", "count_overlap"])
print(mwu, p)

plot.taxon2.unique()
#%%
matched.groupby("arch")["count_overlap"].mean()

#%%
matched.groupby(["mrca_2", "core_remodeling"])["enh_id"].count()
#%%
matched.head()
blankdf = pd.DataFrame({"arch":["simple", "complexenh", "complexenh"],
"count_overlap":[0,0,0], "core_remodeling":[0,1,1], "mrca_2":[0, 0, 0.131],
"taxon2" :[ "Homo sapiens (0)", "Homo sapiens (0)", "Primate (72)"]})
blankdf
