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

#%%
RE ="/dors/capra_lab/projects/enhancer_ages/fantom/results/for_publication/age_breaks/"


colors = [ "amber", "faded green", "dusty purple", "windows blue","greyish"]
palette = sns.xkcd_palette(colors)
sns.palplot(palette)


shuf_colors = [ "amber", "greyish",]
shuf_pal = sns.xkcd_palette(shuf_colors)


#%% Files

PATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/"

enh = "%snon-genic/no-exon_all_fantom_enh/breaks/no-exon_all_fantom_enh_ages_enh_age_arch_summary_matrix.bed" % PATH

shuf = "%snon-genic/no-exon_shuffle_fantom_age_arch_summary_matrix.bed" % PATH
#%% other summary files


# age and taxon file
def get_syn_gen():
    syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
    syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
    syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

    syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2", "mya", "mya2"]] # whittle down the df

    return syn_gen_bkgd


def load_df(f):
    print(f)

    cols = ["#chr_enh", "start_enh", "end_enh", "enh_id", "id",
        "seg_index", "core_remodeling", "arch", "mrca",
        "taxon", "mrca_2", "taxon2", "mya", "mya2"]


    df = pd.read_csv(f, sep = '\t')


    #df = df.rename(columns={'datatype': 'id', 'chr_enh': '#chr_enh'})
    df["enh_len"] = df.end_enh - df.start_enh

    return df


def plot_FigS1_enh_shuf_overlap(shuffle, final_merge):

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

    plt.savefig("%sFigS1.pdf" %RE, bbox_inches = "tight")


def generate_shuf_id(shuffle, final_merge):

    enh_shape = final_merge.shape[0] # get length of enh dataframe

    no_shufs = round(shuffle.shape[0]/enh_shape, 0) # find number of shuffled datasets to label

    sampleD = {} # collect sample results here.

    for n in np.arange(no_shufs):
        sample = shuffle.sample(enh_shape, replace = False) # sample the shuffled dataframe space
        sample["shuf_id"] = n # add a pseudo shuffle id
        sampleD[n] = sample # add sample to the collection dictionary

    new_shuf = pd.concat(sampleD.values())
    return new_shuf


def get_arch_freq(df):

    var = "enh_id"

    if "shuf_id" in list(df):
        cols = ["enh_id", "core_remodeling", "shuf_id"]
        count_cols = ["core_remodeling", "shuf_id"]
        total_col = "shuf_id"

    else:
        cols = ["enh_id", "core_remodeling"]
        count_cols = ["core_remodeling"]

    arch_ = df[cols].drop_duplicates()
    arch_freq = arch_.groupby(count_cols)[var].count().reset_index()

    if "shuf_id" in list(df):
        totals = arch_.groupby([total_col])[var].count().reset_index()
        totals.columns = ["shuf_id", "totals"]
        arch_freq = pd.merge(arch_freq, totals, how = "left")
        arch_freq["dataset"] = "SHUFFLE"

    else:
        totals = arch_.shape[0]
        arch_freq["totals"] = totals
        arch_freq["dataset"] = "FANTOM"

    arch_freq["freq"] = arch_freq[var].divide(arch_freq["totals"])

    return arch_freq


def plot_freq_simple(shuf_arch_freq, arch_freq):

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
    plt.savefig("%ssimple_freq.pdf" %RE, bbox_inches = "tight")


def plot_simple_v_complex(arch_freq):

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

    plt.savefig("%ssimple_v_complex_freq.pdf" %RE, bbox_inches = "tight")


def get_seg_freq(df):

    max_val = "seg_index"
    var = "enh_id"

    if 'shuf_id' in list(df):
        cols = ["enh_id", "shuf_id"]
        freq_cols = ["seg_index", "shuf_id"]

    else:
        cols = ["enh_id"]
        freq_cols = ["seg_index"]

    # quantify the number of ages per  enhancer
    breaks = df.groupby(cols)[max_val].max().reset_index()

    # count enhancers per number of breaks and shuffle ids
    break_freq = breaks.groupby(freq_cols)[var].count().reset_index()

    # get the total number of enhancers
    if 'shuf_id' in list(df):
        totals = break_freq.groupby("shuf_id")[var].sum().reset_index()
        totals.columns = ["shuf_id", "totals"]
        break_freq = pd.merge(break_freq, totals , how = "left", on = "shuf_id")
        break_freq["dataset"] = "Shuffle"

    else:
        totals = df.shape[0]
        break_freq["totals"] = totals
        break_freq["dataset"] = "FANTOM"


    break_freq["freq"] = break_freq[var].divide(break_freq["totals"])

    return break_freq


def get_cumsum(shuf_break_freq):

    cumsum_dict = {}

    for sid in shuf_break_freq.shuf_id.unique():

        testdf = shuf_break_freq.loc[shuf_break_freq.shuf_id == sid]

        cdf= np.cumsum(testdf.freq)/1

        newdf = pd.DataFrame({"cdf": cdf})

        testdf = pd.merge(testdf, newdf, left_index = True, right_index = True)

        cumsum_dict[sid] = testdf

    shufbreak_freq_cdf = pd.concat(cumsum_dict.values())
    shufbreak_freq_cdf.columns=["seg_index", "shuf_id", "enh_id", "totals",
     "dataset", "shuf_freq", "shuf_cdf"]

    return shufbreak_freq_cdf


def plot_cdf(shufbreak_freq_cdf, breaks_freq):

    shuf_cdfplot = shufbreak_freq_cdf.groupby(["seg_index", "dataset"])[["shuf_freq","shuf_cdf", "enh_id", "totals"]].mean().reset_index()

    shuf_cdfplot.columns = ['seg_index', 'dataset', 'freq', 'cdf', "enh_id", "totals"]

    concat = [breaks_freq, shuf_cdfplot]
    plot_cdf_ = pd.concat(concat)


    # the plotting function

    shuf_colors = [ "amber", "greyish",]
    shuf_pal = sns.xkcd_palette(shuf_colors)

    fig, ax = plt.subplots(figsize = (8,8))

    x = "seg_index"
    y = "cdf"
    hue = "dataset"

    sns.lineplot(x = x, y = y, data = plot_cdf_,
    palette = shuf_pal,
    hue = hue,
    )
    ax.set(xticks = (np.arange(0, plot_cdf_.seg_index.max(), step = 1)), \
    xlabel = "number of segments", ylabel = "cumulative distribution",
    xlim = (0,10))

    plt.savefig("%sFantom_CDF_breaks.pdf" %RE, bbox_inches = 'tight')


def get_OR_age_seg(breaks_freq):

    OR_dict = {}

    for seg_index in breaks_freq.seg_index.unique():

        total_enh = breaks_freq.enh_id.sum()
        a = breaks_freq.loc[breaks_freq.seg_index ==seg_index, "enh_id"].tolist()[0] # num simple enhancers
        b = total_enh - a # num complex enhancers
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

    return ORdf


def plot_Fig2B_age_seg_OR(ORdf):
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
# tissue/cell line dataset descriptions
desc_file = "/dors/capra_lab/data/fantom/fantom5/facet_expressed_enhancers/sample_id_descriptions.txt"
desc_df= pd.read_csv(desc_file, sep = '\t', header = None)



#%% LOAD Files


shuffle = load_df(shuf)

final_merge = load_df(enh)

syn_gen_bkgd = get_syn_gen()


#%% are the shuffle and fantom architecture legnths similar?
plot_FigS1_enh_shuf_overlap(shuffle, final_merge)

#%% stratify shuffled datasets from 90x shuffle.

shuffle = generate_shuf_id(shuffle, final_merge)

#%% count the number of simple, complex enhancers in shuffle
# get frequency of simple/ complex enhancers

shuf_arch_freq = get_arch_freq(shuffle)

arch_freq = get_arch_freq(final_merge)

arch_freq
#%% PLOT FANTOM simple v. SHUFFLE simple (64% v. 58% simple enhancers)

plot_freq_simple(shuf_arch_freq, arch_freq)

#%% PLOT FANTOM simple v. COMPLEX (64% simple v. 36% complex enhancers)
plot_simple_v_complex(arch_freq)



#%% get shuffled breaks distribution
shuf_break_freq = get_seg_freq(shuffle)

shuf_break_freq


#%% get shuffled cumulative dist frequency of the breaks.

shufbreak_freq_cdf = get_cumsum(shuf_break_freq)

shufbreak_freq_cdf.head()



#%% get enhancer breaks distribution, CDF
breaks_freq = get_seg_freq(final_merge)

breaks_freq["cdf"]= np.cumsum(breaks_freq.freq)/1
breaks_freq
outf = f"{RE}fantom_num_seg_counts.tsv"
breaks_freq.to_csv(outf, sep = '\t', index = False)


#%% plot cdf

plot_cdf(shufbreak_freq_cdf, breaks_freq)

#%% Are there fewer breaks than expected? Do an FET


archs = pd.merge(shufbreak_freq_cdf, breaks_freq, how = "left", on = "seg_index" )


total_shuf_breaks = shuf_break_freq.groupby(["seg_index"])["enh_id"].sum().reset_index()
total_shuf = total_shuf_breaks.enh_id.sum()
print(total_shuf)

#%%


#%% plot the OR

ORdf = get_OR_age_seg(breaks_freq)

plot_Fig2B_age_seg_OR(ORdf)
