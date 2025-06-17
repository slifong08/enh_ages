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
RE ="/dors/capra_lab/projects/enhancer_ages/roadmap_encode/results/age_breaks/relative_simple/non-genic/"



colors = [ "amber", "faded green", "dusty purple", "windows blue","greyish"]
palette = sns.xkcd_palette(colors)
sns.palplot(palette)

#%%


simple_break_percentile = 0.5

CALCULATE_OR = 1 # if 1, recalculate ORs

#%% Files


pre_path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/"
path = "%sbreaks/non-genic/" % pre_path

# file contains both enhancer + 10 shuffled breaks.
enhFs= glob.glob("%sno-exon_ROADMAP_E*_enh_and_shuf_age_arch_summary_matrix.bed" % (path))

enh_percentiles = "%sbreaks/all_ROADMAP_breaks_percentiles.bed" % pre_path
pdf = pd.read_csv(enh_percentiles, sep = '\t')

len(enhFs)
#%%
# age and taxon file
syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2", "mya", "mya2"]] # whittle down the df
syn_gen_bkgd


# tissue/cell line dataset descriptions
desc_file = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/roadmap_hg19_sample_id_desc.csv"
desc_df= pd.read_csv(desc_file, sep = '\t', header = None)


#%% get relative cut off


def get_relative_simple(sid, pdf):
    if sid in pdf.sid2.unique():
        test = pdf.loc[(pdf.sid2 == sid)]

        if simple_break_percentile in test.pct.unique():
            percentile_used = 0.5
            relative_simple = pdf.loc[(pdf.sid2 == sid) &\
             (pdf.pct == simple_break_percentile), "seg_index"].iloc[0].astype(int)

        else: # idk why some datasets do not have a 0.5 percentile
            relative_simple = pdf.loc[(pdf.sid2 == sid) &\
             (pdf.pct == 0.6), "seg_index"].iloc[0].astype(int)
            percentile_used = 0.6
        return relative_simple, percentile_used
    else:
        return "no_relative_simple", "no_percentile_used"


#%% get odds ratio per break


def get_OR(final_merge, shuffle, seg_index, sid):
    unique_key = sid + "-" + str(seg_index)

    total_enh = final_merge.enh_id.count()
    total_shuf = shuffle.enh_id.count()

    a = final_merge.loc[final_merge.seg_index == seg_index, "enh_id"].count()
    b = total_enh - a
    c = shuffle.loc[shuffle.seg_index == seg_index, "enh_id"].count()
    d = total_shuf - c

    obs = [[a,b], [c,d]]
    OR, P = stats.fisher_exact(obs)
    table = sm.stats.Table2x2(obs) # get confidence interval
    odds_ci = table.oddsratio_confint()

    results = pd.DataFrame({"sid":[sid], "obs": [obs],
    "OR": [OR], "ci":[odds_ci],
    "P": [P], "relative_simple_seg_index": [relative_simple],
     "percentile_used":[percentile_used], "seg_index":[seg_index], "unqkey":[unique_key]})

    return results, unique_key

#%%#%% LOAD Files


if CALCULATE_OR == 1:
    OR_dict = {}
    seg_dict = {}
#%%
if CALCULATE_OR == 1:
    for enh in enhFs:

        sid = (enh.split("/")[-1]).split("_")[2]
        relative_simple, percentile_used = get_relative_simple(sid, pdf)
        print(sid, relative_simple)

        if relative_simple != "no_relative_simple":
            print(sid, relative_simple)
            df = pd.read_csv(enh, sep = '\t', header = None,low_memory=False,
             error_bad_lines = False)

            df = df.loc[df[16]!="datatype"]


            df.columns = ['chr_enh','start_enh','end_enh','shuf_id','enh_id',
            'core_remodeling','arch', 'seg_index', 'mrca','enh_len','taxon','mrca_2',
            'taxon2','mya','mya2','seg_den','datatype']


            # RELATIVE SIMPLE DEF
            df["relative_arch"] = "rel_simple"
            df.loc[df.seg_index.astype(float) >=relative_simple, "relative_arch"] = "rel_complex"
            df.loc[df.relative_arch == "rel_simple", "core_remodeling"] = 0



            # separate dataframe into enhancer and shuffle coordinates.

            if len(df.datatype.unique())>1:
                final_merge = df.loc[~ df.datatype.str.contains("shuf")]

                shuffle = df.loc[df.datatype.str.contains("shuf")]


                # get the frequency of each break in shuffled datasets


                shuf_arch = shuffle[["enh_id", "core_remodeling", "shuf_id"]]
                shuf_arch_freq = shuf_arch.groupby(["core_remodeling", "shuf_id"])["enh_id"].count().reset_index()


                totals = shuf_arch.groupby(["shuf_id"])["enh_id"].count().reset_index()
                totals.columns = ["shuf_id", "totals"]
                shuf_arch_freq = pd.merge(shuf_arch_freq, totals, how = "left")
                shuf_arch_freq["freq"] = shuf_arch_freq["enh_id"].divide(shuf_arch_freq.totals)
                shuf_arch_freq["dataset"] = "SHUFFLE"
                shuf_arch_freq["core_remodeling"] = shuf_arch_freq["core_remodeling"].astype(int)

                # get the frequency of each break in enhancer datasets


                arch = final_merge[["enh_id", "core_remodeling"]].drop_duplicates()
                arch_freq = arch.groupby("core_remodeling")["enh_id"].count().reset_index()
                arch_freq.core_remodeling = arch_freq.core_remodeling.astype(int)
                totals = len(arch.enh_id.unique())

                arch_freq["freq"] = arch_freq["enh_id"].divide(totals)
                arch_freq["dataset"] = "ROADMAP"


                # FET of relative simple v. complex

                # get the total number of simple and complex breaks in the shuffle datasets.
                shuf_sum = shuf_arch_freq.groupby(["core_remodeling"])["enh_id"].sum().reset_index()

                a = arch_freq.loc[arch_freq.core_remodeling ==0, "enh_id"].iloc[0]# num simple enhancers
                b = arch_freq.loc[arch_freq.core_remodeling !=0, "enh_id"].iloc[0]
                c = shuf_sum.loc[shuf_sum.core_remodeling ==0, "enh_id"].iloc[0] # num simple shuffle
                d = shuf_sum.loc[shuf_sum.core_remodeling !=0, "enh_id"].iloc[0] # num simple shuffle

                obs = [[a,b], [c,d]]
                OR, P = stats.fisher_exact(obs)
                table = sm.stats.Table2x2(obs) # get confidence interval
                odds_ci = table.oddsratio_confint()

                results = pd.DataFrame({"sid":[sid], "obs": [obs],
                "OR": [OR], "ci":[odds_ci],
                "P": [P], "relative_simple_seg_index": [relative_simple],
                 "percentile_used":[percentile_used]})

                OR_dict[sid]= results
                #print(obs, OR, P)
                for seg_index in final_merge.seg_index.unique():
                    seg_results, unqkey = get_OR(final_merge, shuffle, seg_index, sid)
                    seg_dict[unqkey] = seg_results
        else:
            print("run shuffles for", sid)


    # save the seg_index, relative_simple OR results


    results = pd.concat(OR_dict.values())

    results["sig"] = 0
    results.loc[results.P<0.05, "sig"] = 1
    results["simple_or"] = 0
    results.loc[results.OR>1, "simple_or"] = 1

    outf = "%sno-exon_ALL_ROADMAP_Shuf_summary_relative_simple_OR.tsv" %path
    results.to_csv(outf, sep ='\t', index = False)


    #  seg_index results


    seg_results = pd.concat(seg_dict.values())

    seg_results["log2OR"] = np.log2(seg_results.OR)

    outf_seg = "%sno-exon_ALL_ROADMAP_Shuf_summary_seg_index_OR.tsv" %path
    seg_results.to_csv(outf_seg, sep ='\t', index = False)

#%%

if CALCULATE_OR ==0:

    outf = "%sno-exon_ALL_ROADMAP_Shuf_summary_relative_simple_OR.tsv" %path
    results = pd.read_csv(outf, sep = '\t')

    outf_seg = "%sno-exon_ALL_ROADMAP_Shuf_summary_seg_index_OR.tsv" %path
    seg_results = pd.read_csv(outf, sep = '\t')


#%%
print(len(results.sid.unique()))
#enhFs.index(enh)

#%%
seg_results.head()


#%% plot per break


seg_results["log2OR"] = np.log2(seg_results.OR)
fig, ax = plt.subplots(figsize =(8,8))

sns.set("poster")
limited_seg_results = seg_results.loc[seg_results.seg_index.astype(float) <11]

sns.barplot(x = "seg_index", y = "log2OR", data = limited_seg_results,
            linewidth=2.5, facecolor=(1, 1, 1, 0),
             edgecolor=".2",)#  yerr=(limited_seg_results["ci"].iloc[0][1] - limited_seg_results["ci"].iloc[0][0]))

ax.set(ylabel= "Fold Change v. Bkgd\n(log2-scaled)",\
 xlabel = "Number of Age Segments",) #xticks = (np.arange(0, limited_seg_results.seg_index.astype(float).max(), step = 10)))#, ylim = (-1.2,0.5))

plt.axhline(0, color = "grey", linewidth = 2.5)

ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(round(2**x, 2)))

ax.yaxis.set_major_formatter(ticks)
ax.yaxis.set_major_locator(MultipleLocator(0.1))

plt.savefig("%sfig2b-ROADMAP_noexon_age_seg_fold_change_matplotlib.pdf" %(RE), bbox_inches = "tight")
#%% how many ORs are significant?


print(results.OR.median())
fig, (ax, ax2) = plt.subplots(ncols = 2, figsize = (12,6))
sns.boxplot(x = "sig", y = "OR", data = results, notch = True, showfliers = True, ax = ax)
sns.swarmplot(x = "sig", y = "OR", data = results, ax = ax)


# how many ORs are significant?
count_df = results.groupby(["sig"])["sid"].count().reset_index()
sigN = count_df.loc[(count_df.sig ==1), "sid"].iloc[0]
not_sigN = count_df.loc[(count_df.sig !=1),  "sid"].iloc[0]
ax_labels = ['not sig. OR\nn=%s' %not_sigN, 'sig. OR\nn=%s'%sigN]

ax.set(title = "all", xticklabels = ax_labels)

# how many ORs are significant and simple enriched?
count_df = results.groupby(["sig", "simple_or"])["sid"].count().reset_index()
simple_sigN = count_df.loc[(count_df.sig ==1)& (count_df.simple_or ==1), "sid"].iloc[0]
complex_sigN = count_df.loc[(count_df.sig ==1)& (count_df.simple_or !=1),  "sid"].iloc[0]
ax2_labels = ['complex\nenriched\nn=%s' % complex_sigN,'simple\nenriched\nn=%s' %simple_sigN]


sns.boxplot(x = "simple_or", y = "OR", data = results.loc[results.sig ==1], ax = ax2)
sns.swarmplot(x = "simple_or", y = "OR", data = results.loc[results.sig ==1], ax = ax2)
ax2.set(title = "sig_only", xticklabels = ax2_labels)
print(results.groupby("simple_or")["sid"].count())
plt.savefig("%sROADMAP_noexon_enh_shuf_simple_OR.pdf" % (RE), bbox_inches = "tight")

print(results.groupby(["sig", "simple_or"])["sid"].count())
#%% PLOT ROADMAP simple v. SHUFFLE simple (64% v. 58% simple enhancers)

archs = pd.concat([shuf_arch_freq, arch_freq]) # combine datasets for plotting

shuf_colors = [ "amber", "greyish",]
shuf_pal = sns.xkcd_palette(shuf_colors)
hue_order = ["ROADMAP", "SHUFFLE"]

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
plt.savefig("%sROADMAP_noexon_enh_shuf_simple_freq.pdf" % (RE), bbox_inches = "tight")


#%% PLOT ROADMAP simple v. COMPLEX (64% simple v. 36% complex enhancers)


fig, ax = plt.subplots(figsize = (8, 8))
sns.set_context("poster")
sns.barplot(x = "core_remodeling", y="freq", data = arch_freq, palette = palette)
ax.set(xticklabels= ["Simple", "Complex\nEnhancer"], xlabel = "", \
ylabel= "Frequency of Dataset", title= "ROADMAP 82/98 Enhancer Architectures")
for p in ax.patches:
    x=p.get_bbox().get_points()[:,0]
    y=p.get_bbox().get_points()[1,1]
    ax.annotate('{:.0f}%'.format(100.*y), (x.mean(), y),
            ha='left', va='bottom', color = "k", alpha = 0.4, fontsize = 20) # set the alignment of the text

plt.savefig("%sROADMAP_noexon_enh_simple_v_complex_freq.pdf" %(RE), bbox_inches = "tight")

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

breaks_freq["dataset"] = "ROADMAP"
breaks_freq["shuf_id"] = "ROADMAP"
breaks_freq["cdf"]= np.cumsum(breaks_freq.freq)/1

breaks_freq.head()


#%% plot cdf


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
xlim = (0,30))
plt.savefig("%sROADMAP_noexon_CDF.pdf" % (RE), bbox_inches = 'tight')

#%%
missing = ["E063", "E105", "E056", "E111", "E012", "E034", "E098", "E029", "E008", "E102", "E078", "E066", "E114", "E103", "E055"]
print(len(missing))
for sid in missing:
    relative_simple, percentile_used = get_relative_simple(sid, pdf)
    print(sid, relative_simple)
