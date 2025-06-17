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


colors = [ "amber", "faded green", "dusty purple", "windows blue","greyish"]
palette = sns.xkcd_palette(colors)
sns.palplot(palette)


shuf_colors = [ "amber", "greyish",]
shuf_pal = sns.xkcd_palette(shuf_colors)


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


#%% are the shuffle and fantom architecture legnths similar?

sns.distplot(shuffle.enh_len.loc[shuffle.arch == "simple"],
color = "grey", kde_kws={'linestyle':'--'},
label = "shuffle-simple")
sns.distplot(shuffle.enh_len.loc[shuffle.arch != "simple"],
color = "k", kde_kws={'linestyle':'--'},
label = "shuffle-complex")
sns.distplot(final_merge.enh_len.loc[final_merge.arch == "simple"],
color = "gold",
 label = "FANTOM-simple")
sns.distplot(final_merge.enh_len.loc[final_merge.arch != "simple"],
color = "g",
label = "FANTOM-complex")
plt.xlabel("length")
plt.ylabel("KDE")
plt.legend()
plt.savefig("%sfantom_enh_shuf_len_dist.pdf" %RE, bbox_inches = "tight")

#%% count the number of simple, complex enhancers in shuffle
shuf_arch = shuffle[["enh_id", "core_remodeling", "shuf_id"]].drop_duplicates()
shuf_arch_freq = shuf_arch.groupby(["core_remodeling", "shuf_id"])["enh_id"].count().reset_index()


# get frequency of simple/ complex enhancers

totals = shuf_arch.groupby(["shuf_id"])["enh_id"].count().reset_index()
totals.columns = ["shuf_id", "totals"]
shuf_arch_freq = pd.merge(shuf_arch_freq, totals, how = "left")
shuf_arch_freq["freq"] = shuf_arch_freq["enh_id"].divide(shuf_arch_freq.totals)
shuf_arch_freq["dataset"] = "SHUFFLE"


# count the number of simple, complex enhancers in shuffle


arch = final_merge[["enh_id", "core_remodeling"]].drop_duplicates()
arch_freq = arch.groupby("core_remodeling")["enh_id"].count().reset_index()
totals = len(arch)
arch_freq["freq"] = arch_freq["enh_id"].divide(totals)
arch_freq["dataset"] = "FANTOM"


#%% PLOT FANTOM simple v. SHUFFLE simple (64% v. 58% simple enhancers)


archs = pd.concat([shuf_arch_freq, arch_freq]) # combine datasets for plotting

hue_order = ["FANTOM", "SHUFFLE"]

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
plt.savefig("%sfantom_enh_shuf_simple_freq.pdf" %RE, bbox_inches = "tight")



#%% PLOT FANTOM simple v. COMPLEX (64% simple v. 36% complex enhancers)



fig, ax = plt.subplots(figsize = (8, 8))
sns.set_context("poster")
sns.barplot(x = "core_remodeling", y="freq", data = arch_freq, palette = palette)
ax.set(xticklabels= ["Simple", "Complex\nEnhancer"], xlabel = "", \
ylabel= "Frequency of Dataset", title= "FANTOM Enhancer Architectures")
for p in ax.patches:
    x=p.get_bbox().get_points()[:,0]
    y=p.get_bbox().get_points()[1,1]
    ax.annotate('{:.0f}%'.format(100.*y), (x.mean(), y),
            ha='left', va='bottom', color = "k", alpha = 0.4, fontsize = 20) # set the alignment of the text

plt.savefig("%sfantom_enh_simple_v_complex_freq.pdf" %RE, bbox_inches = "tight")



#%% get shuffled breaks distribution


# quantify the number of ages per shuffled enhancer
shuf_break = shuffle.groupby(["enh_id", "shuf_id"])["seg_index"].max().reset_index()

# count enhancers per number of breaks and shuffle ids
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


#%% get shuffled cumulative dist frequency of the breaks.


cumsum_dict = {}
for sid in shuf_break_freq.shuf_id.unique():
    testdf = shuf_break_freq.loc[shuf_break_freq.shuf_id == sid]
    newdf = get_cumsum(testdf)
    cumsum_dict[sid] = newdf

shufbreak_freq_cdf = pd.concat(cumsum_dict.values())

shufbreak_freq_cdf.columns=["seg_index", "shuf_id", "shuf_count", "shuf_totals", "shuf_freq", "shuf_dataset", "shuf_cdf"]
shufbreak_freq_cdf.head()



#%% get enhancer breaks distribution, CDF

breaks = final_merge.groupby(["enh_id"])["seg_index"].max().reset_index()
breaks_freq = breaks.groupby("seg_index")["enh_id"].count().reset_index()
total_enh = len(breaks)
print(total_enh)

breaks_freq["freq"] = breaks_freq["enh_id"].divide(total_enh)

breaks_freq["dataset"] = "FANTOM"
breaks_freq["shuf_id"] = "FANTOM"
breaks_freq["cdf"]= np.cumsum(breaks_freq.freq)/1

breaks_freq.head()


#%% plot cdf


shuf_cdfplot = shufbreak_freq_cdf[["seg_index", "shuf_freq", "shuf_dataset", "shuf_id", "shuf_cdf"]]
shuf_cdfplot.columns = ['seg_index', 'freq', 'dataset', 'shuf_id', 'cdf']



concat = [breaks_freq, shuf_cdfplot]
plot_cdf = pd.concat(concat)


# the plotting function


fig, ax = plt.subplots(figsize = (8,8))
x = "seg_index"
y = "cdf"

sns.lineplot(x, y, data = plot_cdf, ax = ax, hue = "dataset", palette = shuf_pal)
ax.set(xticks = (np.arange(0, plot_cdf.seg_index.max(), step = 5)), \
xlabel = "number of segments", ylabel = "cumulative distribution",
xlim = (0,10))
plt.savefig("%sFantom_CDF_breaks.pdf" %RE, bbox_inches = 'tight')


#%% Are there fewer breaks than expected? Do an FET


archs = pd.merge(shufbreak_freq_cdf, breaks_freq, how = "left", on = "seg_index" )


total_shuf_breaks = shuf_break_freq.groupby(["seg_index"])["enh_id"].sum().reset_index()
total_shuf = total_shuf_breaks.enh_id.sum()
print(total_shuf)

#%%


OR_dict = {}

for seg_index in breaks_freq.seg_index.unique():
    a = breaks_freq.loc[breaks_freq.seg_index ==seg_index, "enh_id"].tolist()[0] # num simple enhancers
    b = total_enh - a # num complex enhancers
    #b = breaks_freq.loc[breaks_freq.seg_index !=seg_index, "enh_id"].sum() # num simple enhancers
    c = total_shuf_breaks.loc[total_shuf_breaks.seg_index ==seg_index, "enh_id"].tolist()
    if len(c)>0:
        c = c[0]
    else:
        c = 0 # num simple shuffle # num simple shuffle
    d = total_shuf - c # num complex shuffle


    obs = [[a,b], [c,d]]
    OR, P = stats.fisher_exact(obs)
    table = sm.stats.Table2x2(obs) # get confidence interval
    odds_ci = table.oddsratio_confint()
    newdf = pd.DataFrame({"seg_index":[seg_index], "a":[a], "b":[b], "c":[c], "d":[d],
                         "OR":[OR], "P":[P], "ci_lower" :[odds_ci[0]],
                        "ci_upper" :[odds_ci[1]]})
    OR_dict[seg_index] = newdf
    print(seg_index, obs, OR, P)


ORdf = pd.concat(OR_dict.values())
ORdf["yerr"] = ORdf.ci_upper-ORdf.ci_lower
ORdf['log'] = np.log2(ORdf.OR)
ORdf.head()
#%%
a, b
#%%
breaks_freq.head()

#%% plot the OR


fig, ax = plt.subplots(figsize =(8,8))
sns.set("poster")

firstfive = ORdf.loc[ORdf.seg_index <6]

sns.barplot(x = "seg_index", y = "log", data = firstfive.loc[firstfive.seg_index <6],
            linewidth=2.5, facecolor=(1, 1, 1, 0),
             edgecolor=".2",  yerr=(firstfive["ci_upper"] - firstfive["ci_lower"]))

ax.set(ylabel= "Fold Change v. Bkgd\n(log2-scaled)",\
 xlabel = "Number of Age Segments", ylim = (-1.2,0.5))

plt.axhline(0, color = "grey", linewidth = 2.5)

ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(round(2**x, 1)))

ax.yaxis.set_major_formatter(ticks)
ax.yaxis.set_major_locator(MultipleLocator(0.5))


plt.savefig("%sfig2b-fantom_age_seg_fold_change_matplotlib.pdf" %RE, bbox_inches = "tight")

#%%
hue_order = ["FANTOM", "Shuffle"]
fig, ax = plt.subplots(figsize = (8,8))
sns.set("poster")
sns.barplot(x = "seg_index", y = "cdf", data = plot_cdf,
hue = "dataset", hue_order = hue_order, palette = shuf_pal)
ax.set(xlabel="Number of Age Segments per Enhancer",\
ylabel = "Cumulative Frequency of Dataset",\
title = "FANTOM Enhancer Age Segment Count\n Cumulative distribution",
xlim=(-1, 5.5))

ax.set_title("")
ax.get_legend().remove()

plt.savefig("%sfig2b-fantom_age_seg_cum_dist_matplotlib.pdf" %RE, bbox_inches = "tight")
