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
RE ="/dors/capra_lab/projects/enhancer_ages/emera16/results/age_breaks/"
colors = [ "amber", "faded green", "dusty purple", "windows blue","greyish"]
palette = sns.xkcd_palette(colors)
sns.palplot(palette)


#%% Files

path = "/dors/capra_lab/projects/enhancer_ages/emera16/data/non-genic/"

enh = "%sno-exon_emera_2016_neocortical_dev_enhancers_hu_ms_parallel_breaks.bed" % path

shufs = glob.glob("%sno-exon_shuf-emera_2016_neocortical_dev_enhancers_hu_ms-*_parallel_breaks.bed" % path)

#%% other summary files

# age and taxon file
syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2", "mya", "mya2"]] # whittle down the df


#%% LOAD Files
cols = ["chr_enh", "start_enh", "end_enh", "shuf_id",
 "enh_id", "seg_index", "core_remodeling", "arch", "mrca_2"]
shuf_dict = {}
for shuf in shufs:
    shuf_id = (shuf.split("/")[-1]).split(".")[0]

    shuffle = pd.read_csv(shuf, sep = '\t', header =None)
    shuffle.columns = cols
    shuffle["shuf_id"] =shuf_id

    shuf_dict[shuf_id] = shuffle

shuffle = pd.concat(shuf_dict.values())
shuffle = shuffle.loc[shuffle.mrca_2 != "max_age"]
shuffle["mrca_2"] = shuffle["mrca_2"].astype(float).round(3) # round ages

final_merge = pd.read_csv(enh, sep = '\t', header =None)
final_merge.columns = cols
final_merge.shuf_id = "Emera16"
final_merge = final_merge.loc[final_merge.mrca_2 != "max_age"]
final_merge["mrca_2"] = final_merge["mrca_2"].astype(float).round(3) # round ages

relative_simple = final_merge.seg_index.median()

final_merge["core_remodeling"]= 0
shuffle["core_remodeling"]= 0
final_merge.seg_index=final_merge.seg_index.astype(int)
shuffle.seg_index=shuffle.seg_index.astype(int)
final_merge.loc[final_merge.seg_index.astype(int)>=relative_simple, "core_remodeling"] = 1
shuffle.loc[shuffle.seg_index.astype(int) >=relative_simple, "core_remodeling"] = 1
print(relative_simple)
#%%

#%% LOAD Files


shuf_arch = shuffle[["enh_id", "core_remodeling", "shuf_id"]].drop_duplicates()
shuf_arch_freq = shuf_arch.groupby(["core_remodeling", "shuf_id"])["enh_id"].count().reset_index()
shuf_arch_freq.head()

totals = shuf_arch.groupby(["shuf_id"])["enh_id"].count().reset_index()
totals.columns = ["shuf_id", "totals"]
shuf_arch_freq = pd.merge(shuf_arch_freq, totals, how = "left")
shuf_arch_freq["freq"] = shuf_arch_freq["enh_id"].divide(shuf_arch_freq.totals)
shuf_arch_freq["dataset"] = "SHUFFLE"

arch = final_merge[["enh_id", "core_remodeling"]].drop_duplicates()
arch_freq = arch.groupby("core_remodeling")["enh_id"].count().reset_index()
totals = len(arch)
arch_freq["freq"] = arch_freq["enh_id"].divide(totals)
arch_freq["dataset"] = "emera16"


#%% PLOT emera16 simple v. SHUFFLE simple (64% v. 58% simple enhancers)
archs = pd.concat([shuf_arch_freq, arch_freq]) # combine datasets for plotting

shuf_colors = [ "amber", "greyish",]
shuf_pal = sns.xkcd_palette(shuf_colors)
hue_order = ["emera16", "SHUFFLE"]

fig, ax = plt.subplots(figsize = (8, 8))
sns.set_context("poster")
sns.barplot(x = "core_remodeling", y="freq", data = archs.loc[archs.core_remodeling ==0],
            hue = "dataset", hue_order = hue_order, palette = shuf_pal)

for p in ax.patches:
    x=p.get_bbox().get_points()[:,0]
    y=p.get_bbox().get_points()[1,1]
    ax.annotate('{:.0f}%'.format(100.*y), (x.mean(), y),
            ha='left', va='bottom', color = "k", alpha = 0.4, fontsize = 20) # set the alignment of the text

ax.set(xticklabels = "", xlabel= "", title= "", ylabel= "Frequency of Dataset")
ax.get_legend().remove()
plt.savefig("%ssfig2.14_emera16_enh_shuf_simple_freq_rel_simple.pdf" %RE, bbox_inches = "tight")

#%% PLOT emera16 simple v. COMPLEX (64% simple v. 36% complex enhancers)

fig, ax = plt.subplots(figsize = (8, 8))
sns.set_context("poster")
sns.barplot(x = "core_remodeling", y="freq", data = arch_freq, palette = palette)
ax.set(xticklabels= ["Simple", "Complex\nEnhancer"], xlabel = "", \
ylabel= "Frequency of Dataset", title= "emera16 Enhancer Architectures")
for p in ax.patches:
    x=p.get_bbox().get_points()[:,0]
    y=p.get_bbox().get_points()[1,1]
    ax.annotate('{:.0f}%'.format(100.*y), (x.mean(), y),
            ha='left', va='bottom', color = "k", alpha = 0.4, fontsize = 20) # set the alignment of the text

plt.savefig("%ssfig2.14_emera16_enh_simple_v_complex_freq_rel_simple.pdf" %RE, bbox_inches = "tight")

#%% get shuffled breaks distribution
shuf_break = shuffle.groupby(["enh_id", "shuf_id"])["seg_index"].max().reset_index()

shuf_break_freq = shuf_break.groupby(["seg_index", "shuf_id"])["enh_id"].count().reset_index()

shuf_totals = shuf_break_freq.groupby("shuf_id")["enh_id"].sum().reset_index()
shuf_totals.columns = ["shuf_id", "shuf_id_totals"]

shuf_break_freq = pd.merge(shuf_break_freq, shuf_totals, how = "left", on = "shuf_id")


shuf_break_freq["freq"] = shuf_break_freq["enh_id"].divide(shuf_break_freq.shuf_id_totals)
shuf_break_freq["dataset"] = "Shuffle"

#%% cumsum FUNCTION
def get_cumsum(df):

    cdf= np.cumsum(df.freq)/1
    newdf = pd.DataFrame({"cdf": cdf})
    testdf = pd.merge(df, newdf, left_index = True, right_index = True)
    return testdf

#%% get shuffled cumulative dist
cumsum_dict = {}
for sid in shuf_break_freq.shuf_id.unique():
    testdf = shuf_break_freq.loc[shuf_break_freq.shuf_id == sid]
    newdf = get_cumsum(testdf)
    cumsum_dict[sid] = newdf

shufbreak_freq_cdf = pd.concat(cumsum_dict.values())

shufbreak_freq_cdf.columns=["seg_index", "shuf_id", "shuf_count", "shuf_totals", "shuf_freq", "shuf_dataset", "shuf_cdf"]
shufbreak_freq_cdf.head()

#%% get enhancer breaks distribution
breaks = final_merge.groupby(["enh_id"])["seg_index"].max().reset_index()
breaks_freq = breaks.groupby("seg_index")["enh_id"].count().reset_index()
enh_totals = len(breaks)

breaks_freq["freq"] = breaks_freq["enh_id"].divide(enh_totals)

breaks_freq["dataset"] = "emera16"
breaks_freq["shuf_id"] = "emera16"
breaks_freq["cdf"]= np.cumsum(breaks_freq.freq)/1
#breaks_freq = breaks_freq.drop(["enh_id"], axis = 1)
breaks_freq.head()
#%%
list(breaks_freq)
#%% make shuffle cdf look like breaks_freq
shuf_cdfplot = shufbreak_freq_cdf[["seg_index", "shuf_freq", "shuf_dataset", "shuf_id", "shuf_cdf"]]
shuf_cdfplot.columns = ['seg_index', 'freq', 'dataset', 'shuf_id', 'cdf']
#%%
concat = [breaks_freq, shuf_cdfplot]
plot_cdf = pd.concat(concat)
fig, ax = plt.subplots(figsize = (8,8))
x = "seg_index"
y = "cdf"


sns.lineplot(x, y, data = plot_cdf, ax = ax, hue = "dataset", palette = shuf_pal)
ax.set(xticks = (np.arange(0, plot_cdf.seg_index.max(), step = 5)), \
xlabel = "number of segments", ylabel = "cumulative distribution",
xlim = (0,31))
ax.axvline(relative_simple)
plt.savefig("%ssfig2.14_emera16_CDF_breaks_rel_simple.pdf" %RE, bbox_inches = 'tight')

#%% Are there fewer breaks than expected? Do an FET
archs = pd.merge(shufbreak_freq_cdf, breaks_freq, how = "left", on = "seg_index" )
archs.seg_index = archs.seg_index.astype(int)

total_shuf_breaks = shuf_break_freq.groupby(["seg_index"])["enh_id"].sum().reset_index()

total_shuf_breaks.loc[total_shuf_breaks.seg_index ==1, "enh_id"][0]


#%%
OR_dict = {}

for seg_index in breaks_freq.seg_index.unique():

    a = breaks_freq.loc[breaks_freq.seg_index ==seg_index, "enh_id"].to_list()[0] # num simple enhancers
    b = enh_totals - a # num complex enhancers
    c = total_shuf_breaks.loc[total_shuf_breaks.seg_index ==seg_index, "enh_id"].to_list()
    if len(c)>0:
        c = c[0]
    else:
        c = 0 # num simple shuffle
    d = total_shuf_breaks.enh_id.sum() - c # num complex shuffle

    obs = [[a,b], [c,d]]
    OR, P = stats.fisher_exact(obs)
    table = sm.stats.Table2x2(obs) # get confidence interval
    odds_ci = table.oddsratio_confint()
    newdf = pd.DataFrame({"seg_index":[seg_index], "a":[a], "b":[b], "c":[c], "d":[d],
                         "OR":[OR], "P":[P], "ci_lower" :[odds_ci[0]],
                        "ci_upper" :[odds_ci[1]]})
    OR_dict[seg_index] = newdf
    #print(seg_index, obs, OR, P)

ORdf = pd.concat(OR_dict.values())
ORdf["yerr"] = ORdf.ci_upper-ORdf.ci_lower
ORdf['log'] = np.log2(ORdf.OR)
ORdf.head()
#%%
fig, ax = plt.subplots(figsize =(16,8))
sns.set("poster")
max_segs =25
ORdf.seg_index = ORdf.seg_index.astype(int)
firstfive = ORdf.loc[ORdf.seg_index.astype(int) <max_segs]

sns.barplot(x = "seg_index", y = "log", data = firstfive.loc[firstfive.seg_index <max_segs],
            linewidth=2.5, facecolor=(1, 1, 1, 0),
             edgecolor=".2",  yerr=(firstfive["ci_upper"] - firstfive["ci_lower"]))

ax.set(ylabel= "Fold Change v. Bkgd\n(log2-scaled)",\
 xlabel = "Number of Age Segments")
 #, ylim = (-1.2,0.5))

plt.axhline(0, color = "grey", linewidth = 2.5)

ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(round(2**x, 1)))

ax.yaxis.set_major_formatter(ticks)
ax.yaxis.set_major_locator(MultipleLocator(0.5))



plt.savefig("%sfig2.14-emera16_age_seg_fold_change_matplotlib_rel_simple.pdf" %RE, bbox_inches = "tight")

#%%
plot_cdf.head()
#%%
final_merge.head()
shuffle.head()
#%%
def get_arch_freq(df):

    age_cols = ["core_remodeling", "mrca_2", "mrca_count"]
    age_counts = df.groupby(["core_remodeling", "mrca_2"])["enh_id"].count().reset_index()
    age_counts.columns= age_cols


    arch_cols = ["core_remodeling", "arch_total_count"]
    arch_totals = age_counts.groupby("core_remodeling")["mrca_count"].sum().reset_index()
    arch_totals.columns = arch_cols

    age_freq = pd.merge(age_counts, arch_totals, how = "left", on = "core_remodeling")

    age_freq["arch_freq"] = age_freq["mrca_count"].divide(age_freq["arch_total_count"])

    return age_freq
#%%
enhAgeFreq = get_arch_freq(final_merge)
shufAgeFreq = get_arch_freq(shuffle)
enhAgeFreq.head()

#%%

Xticklabels = ["Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Vert"]
fig, (ax, ax2) = plt.subplots(ncols = 2, figsize = (16,6))
sns.barplot(x = "mrca_2", y = 'arch_freq', data = enhAgeFreq,
hue = "core_remodeling", palette = palette, ax = ax)

ax.set(xticklabels =Xticklabels, ylabel = '% of architecture',
title = "Emera 16 architecture frequency", xlabel = "")

ax.legend().remove()

ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)

sns.barplot(x = "mrca_2", y = 'mrca_count', data = enhAgeFreq,
hue = "core_remodeling", palette = palette, ax = ax2)

ax2.set(xticklabels =Xticklabels, ylabel = 'count',
title = "Emera 16 architecture count", xlabel = "")

ax2.legend().remove()

ax2.set_xticklabels(ax2.get_xticklabels(), rotation = 90)
plt.tight_layout()
plt.savefig("%sfig2.14-emera16_mrca_x_arch_rel_simple.pdf" %RE, bbox_inches = "tight")

#%%
shufcols = ["core_remodeling", "mrca_2", "shuf_mrca_count", "shuf_arch_total", "shuf_arch_freq"]
shufAgeFreq.columns = shufcols
ratio = pd.merge(enhAgeFreq, shufAgeFreq, how = "left", on = ["core_remodeling", "mrca_2"])
ratio["obs_over_exp"] = ratio["arch_freq"].divide(ratio["shuf_arch_freq"])
ratio["log2oe"] = np.log2(ratio["obs_over_exp"])
ratio.head()
#%%
Xticklabels = ["Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Vert"]
fig, (ax) = plt.subplots(ncols = 1, figsize = (6,6))
sns.barplot(x = "mrca_2", y = 'log2oe', data = ratio,
hue = "core_remodeling", palette = palette, ax = ax)

ax.set(xticklabels =Xticklabels, ylabel = 'Fold change\n(log2-scaled)',
title = "Emera 16 architecture fold-change", xlabel = "")

ax.legend().remove()

ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)

plt.tight_layout()
plt.savefig("%sfig2.14-emera16_mrca_fold_change_rel_simple.pdf" %RE, bbox_inches = "tight")
