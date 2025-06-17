
from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy import stats
import seaborn as sns
import statsmodels
import statsmodels.api as sm
import subprocess

FANTOMPATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/non-genic/ages/"
FANTOMFILE = "syn_breaks_no-exon_all_fantom_enh_ages.bed"
FANTOM = os.path.join(FANTOMPATH, FANTOMFILE)

ENCODEPATH = "/dors/capra_lab/data/encode/encode3_hg38/TF/liftOver_hg19"
ENCODEFILE = "trimmed_encRegTfbsClusteredWithCells.liftOver.to.hg19.bed"
ENCODE = os.path.join(ENCODEPATH, ENCODEFILE)

INTERSECTIONPATH = FANTOMPATH
INTERSECTIONFILE = "All_FANTOM_x_ENCODE_hg19.bed"
INTERSECTION = os.path.join(INTERSECTIONPATH, INTERSECTIONFILE)

RE = "/dors/capra_lab/projects/enhancer_ages/fantom/results/encode3/"
RE = "/dors/capra_lab/projects/enhancer_ages/landscape/results/fantom/encode3/"
#%% Functions

def bed_intersect(fantom, encode, intersection):

    if os.path.exists(intersection) == False:
        cmd = "bedtools intersect -a %s -b %s -wao > %s" % (fantom, encode, intersection)

        subprocess.call(cmd, shell = True)

        print(cmd)
    else:
        print("previously done enh x encode intersection")


def load_syn_gen_bkgd(build):

    F = f"/dors/capra_lab/projects/enhancer_ages/{build}_syn_taxon.bed"
    syngenbkgd = pd.read_csv(F, sep='\t')
    syngenbkgd[["mrca", "mrca_2"]] = syngenbkgd[["mrca", "mrca_2"]].round(3)

    return syngenbkgd


def get_core_age(df):
    core = df.groupby("enh_id")["mrca_2"].max().reset_index()
    core.columns = ["enh_id", 'core_mrca_2']
    df = pd.merge(df, core, how = "left")

    return df


def format_df(df):

    # add architecture label
    df["arch"] = "simple"
    df.loc[(df.core_remodeling ==1) & (df.core ==1), "arch"] = "complex_core"
    df.loc[(df.core_remodeling ==1) & (df.core ==0), "arch"] = "complex_derived"

    df["core_remodeling_2"] = 0
    df.loc[df["arch"] == "complex_core", "core_remodeling_2"] = 1
    df.loc[df["arch"] == "complex_derived", "core_remodeling_2"] = 2

    df["overallarch"] = "simple"
    df.loc[df.core_remodeling ==1, "overallarch"] = "complex"

    # add syn identifier
    df["syn_id"] = df.chr_syn + ":" + df.start_syn.map(str) + "-" + df.end_syn.map(str)
    df["tf"] = df.tf_id.apply(lambda x: x.split("_")[0])
    #calculate enhancer and syntenic block length.

    df["enh_len"] = df.end - df.start
    df["syn_len"] = df.end_syn - df.start_syn

    # binary for TF overlap
    df["tfoverlap_bin"] = 1
    df.loc[df.tf == ".", "tfoverlap_bin"] = 0
    df.loc[df.overlap <6, "tfoverlap_bin"] = 0

    build = "hg19"
    syn_gen_bkgd = load_syn_gen_bkgd(build)

    df["mrca"] = df["mrca"].round(3)
    df = pd.merge(df, syn_gen_bkgd, how = "left")
    # get the core age
    df = get_core_age(df)
    df["counts"] = 1
    return df


def make_pdf(file_name, RE):

    OUTFILE = file_name + ".pdf"
    OUTF = os.path.join(RE, OUTFILE)

    return OUTF


def prep_2x2(tf, arch, df):


    dfarch = df.loc[df.arch == arch]
    dfbkgd = df.loc[df.arch != arch]


    TF_in_arch = len(dfarch.loc[dfarch.tf == tf])
    TF_bkgd = len(dfbkgd.loc[dfbkgd.tf == tf])
    not_TF_in_arch = len(dfarch.loc[dfarch.tf != tf])
    not_TF_in_bkgd = len(dfbkgd.loc[dfbkgd.tf != tf])

    a, b, c, d = TF_in_arch, not_TF_in_arch, TF_bkgd, not_TF_in_bkgd

    obs = [[a,b], [c,d]]

    return obs


def quantify_2x2(obs, comparison_name):

    if obs[0][0] > 100:

        OR, P = stats.fisher_exact(obs)
        table = sm.stats.Table2x2(obs) # get confidence interval
        odds_ci = table.oddsratio_confint()
        newdf = pd.DataFrame({"comparison_name":comparison_name,
                              "a":obs[0][0], "b":obs[0][1],
                              "c":obs[1][0], "d":obs[1][1],
                              "OR":[OR], "P":[P],
                              "ci_lower" :[odds_ci[0]],
                              "ci_upper" :[odds_ci[1]],
                            })

        if P<0.05:
            print(comparison_name, obs, round(OR, 2), round(P, 4))

        return newdf


def fdr_correction(collection_dict):

    df = pd.concat(collection_dict.values())

    pvals = df["P"]

    df["reject_null"], df["FDR_P"] = statsmodels.stats.multitest.fdrcorrection(pvals, alpha=0.05)

    return df

def get_tfbs_overlap_fraction(df, groupby, val):

    overlap = df.loc[df.tfoverlap_bin > 0].groupby(groupby)[val].count().reset_index()
    overlap.columns = ["arch", "overlap"]

    nooverlap = df.loc[df.tfoverlap_bin == 0].groupby(groupby)[val].count().reset_index()
    nooverlap.columns = ["arch", "nooverlap"]

    overlap = pd.merge(overlap, nooverlap, how  = "left")

    overlap["total"] = overlap["overlap"] + overlap["nooverlap"]
    overlap["frac_overlap"] = overlap.overlap.divide(overlap.total)
    overlap["frac_nooverlap"] = overlap.nooverlap.divide(overlap.total)

    outfile = f"{RE}fraction_TFBS_overlap_{groupby}.tsv"
    overlap.to_csv(outfile, sep = '\t')

    return overlap


def make_empty_dict(df, groupby_cols):
    # for when i need an empty dictionary with a complete set of architectures and values
    emptydf_dict = {}
    val = 0 # unique identifier
    if len(groupby_cols) == 2:
        for mrca_2 in df[groupby_cols[0]].unique():
            for arch in df[groupby_cols[1]].unique():
                emptydf = pd.DataFrame({ groupby_cols[0]:[mrca_2], groupby_cols[1]:[arch],})
                emptydf_dict[val] = emptydf
                val+=1

    elif len(groupby_cols) == 1:
        for mrca_2 in df[groupby_cols[0]].unique():
            emptydf = pd.DataFrame({ groupby_cols[0]:[mrca_2]})
            emptydf_dict[val] = emptydf
            val+=1
    empty = pd.concat(emptydf_dict.values())

    return empty


def get_counts(df, groupby_cols, groupby_val):

    counts = df.groupby(groupby_cols)[groupby_val].sum().reset_index()


    empty = make_empty_dict(df, groupby_cols) # make an empty df to fill architectures w/ no counts
    counts = pd.merge(empty, counts, how = "left").fillna(0)

    # change count data to int
    counts[groupby_val] = counts[groupby_val].astype(int)

    # sort and reset the index. Seaborn plots by index value.
    counts = counts.sort_values(by = groupby_cols).reset_index()

    # drop the index column.
    counts = counts.drop(["index"], axis = 1)

    return counts


def plot_annotate_counts(splot, counts_df, groupby_val, height_adjust):

    # annotate plot with counts
    for n, p in enumerate(splot.patches):

        value = counts_df.iloc[n][groupby_val].astype(int)
        #print(n, p, value)
        if height_adjust == 0:
            height_adjust = (p.get_height()-0.03)


        splot.annotate(value,
                       (p.get_x() + p.get_width() / 2.,height_adjust),
                       ha = 'center', va = 'baseline',
                       size=15,
                       rotation = 90,
                       color = "k",
                       xytext = (0, 1),
                       textcoords = 'offset points'
                       )

def plot_figure3(x, y, df, fig_id, re, trim_len, frac, dataset):
    # sns colors
    arch_colors = [ "amber", "dusty purple", "windows blue","greyish"]
    arch_palette = sns.xkcd_palette(arch_colors)

    pal = [ "amber", "faded green"]
    palette = sns.xkcd_palette(pal)

    # for ranking and stratifying by age and architecture, set x1 for counts
    if x == 'overallarch':
        x1 = "core_remodeling"
        pal = palette
        order_labs = ["simple", "complex"]
        order = [0,1]
    elif x == "arch":
        x1 = "core_remodeling_2"
        pal = arch_palette
        order_labs = ["simple", "core", "derived"]
        order = [0,1,2]


    ylab = y

    xlab = ['Homo', 'Prim', 'Euar', 'Bore', 'Euth', 'Ther', 'Mam',
    'Amni', 'Tetr', 'Vert']

    title = dataset

    # set up plot
    sns.set("poster")
    fig = plt.figure(figsize = (12, 8))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 3])
    ax = plt.subplot(gs[0])

    # plot the first panel
    data = df

    splot = sns.barplot(x1, y, data = data,
                palette = pal, order = order,
                ax = ax)

    #ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.set_xticklabels(order_labs, rotation = 90)

    # get counts for annotation
    groupby_cols, groupby_val = [x1], "counts"
    countdf = get_counts(df, groupby_cols, groupby_val)
    height_adjust = 0.001
    plot_annotate_counts(splot, countdf, groupby_val, height_adjust)

    # plot the second panel
    ax2 = plt.subplot(gs[1])

    x = "core_mrca_2"
    data = df.sort_values(by = "core_mrca_2")
    hue = x1

    # plot mrca-stratified barplot

    mplot = sns.barplot(x, y, data = data,
                hue = hue,
                palette = pal,
                ax = ax2)

    # add counts and annotate barplot
    groupby_cols, groupby_val = [x1, "core_mrca_2"], "counts"
    countdf = get_counts(df, groupby_cols, groupby_val)
    plot_annotate_counts(mplot, countdf, groupby_val, height_adjust)

    # plot x labels
    ax2.set_xticklabels(xlab, rotation = 90)
    ax2.set(xlabel = "", ylabel = "", title = title)

    ax2.legend().remove()

    #ax2.yaxis.set_major_locator(MultipleLocator(0.1))

    ax2lim = ax2.get_ylim()
    # set first plot lim based on second plot lim
    ax.set(xlabel = "", ylabel = ylab, ylim = ax2lim)

    outf = f"{re}fig{fig_id}_FANTOM_X_ENCODE_{trim_len}_noexon_{dataset}_{frac}_mrcas.pdf"

    plt.savefig(outf, bbox_inches = "tight")


    # do stats
    if x == "arch":
        simple = df.loc[df[x].str.contains('simple'), y]
        complexcore = df.loc[df[x].str.contains('complexcore'), y]
        derived = df.loc[df[x].str.contains('derived'), y]
        k, kp = stats.kruskal(simple,complexcore,derived)
        m, mp = stats.mannwhitneyu(complexcore,  derived)

        stat_df = pd.DataFrame({
        "test": ["krusal-wallis", "MWU"],
        "comparison":["simple v.complexcore v. complexderived", "complexcore v. derived"],
        "stat":[k,m],
        "p-value":[kp, mp]
        })
        medians = df.groupby(x)[y].median().reset_index()
        medians.columns = [x, f"median_{y}"]
        means =  df.groupby(x)[y].mean().reset_index()
        means.columns = [x, f"mean_{y}"]

        out_stat_df = f"{RE}{dataset}_{y}_median_mean.tsv"
        metrics = pd.concat([medians, means])
        metrics.to_csv(out_stat_df, sep = '\t', index = False) # write the stats results

        out_stat_df = f"{RE}{dataset}_{y}_mwu.tsv"
        stat_df.to_csv(out_stat_df, sep = '\t', index = False) # write the stats results



#%% Bed command

bed_intersect(FANTOM, ENCODE, INTERSECTION)

#%% dataframe

cols = ["chr_syn", "start_syn", "end_syn",
"enh_id","chr", "start", "end",
"seg_index", "core_remodeling", "core",
"mrca",
"chr_tf", "start_tf", "end_tf",
"tf_id", "peak_len", "tf",
"cell_line",  "overlap"
]

df_ = pd.read_csv(INTERSECTION,
sep = '\t',
)

df_.columns = cols # add column names


df = format_df(df_) # format the dataframe

df.info()

df.head()



#%% get some basic info about Fantom enhancer overlap


enh_df = df.groupby(["enh_id", "core_remodeling", "overallarch"])[["mrca", "seg_index"]].max().reset_index()

totalenh_n = len(enh_df) #30279 enhancers total
simpleenh_n = len(enh_df.loc[enh_df.overallarch == "simple"]) #14098 simple enhancers
complexenh_n = len(enh_df.loc[enh_df.overallarch != "simple"]) # 8744 complex enhancers


#%% calculate enhancer TF density


tf_density = df.groupby(["enh_id", "enh_len", "overallarch"])["tfoverlap_bin"].sum().reset_index().drop_duplicates()

tf_density["tf_density"] = tf_density["tfoverlap_bin"].divide(tf_density.enh_len)


tf_density = pd.merge(tf_density, df[["enh_id", "core_remodeling", "core_mrca_2", "counts"]], how = "left").drop_duplicates()
tf_density.shape
#%% how many enhancers do not overlap TFs?

groupby = "overallarch"
val = "enh_id"
overlap = get_tfbs_overlap_fraction(tf_density, groupby, val)
overlap
#%%
x, y = "overallarch", "tf_density"

fig_id = "4A"
re = RE
trim_len = "raw"
frac = "all"
dataset = "all_fantom_enh"

tf_density.head()
plot_figure3(x, y, tf_density, fig_id, re, trim_len, frac, dataset)

#%%
x, y = "overallarch", "tf_density"

fig_id = "4A"
re = RE
trim_len = "raw"
frac = "all"
dataset = "all_fantom_enh_non_zero"


plot_figure3(x, y, tf_density.loc[tf_density.tf_density>0], fig_id, re, trim_len, frac, dataset)

#%%


syn_tf_density = df.groupby(["syn_id", "syn_len", "arch"])["tfoverlap_bin"].sum().reset_index()
syn_tf_density["tf_density"] = syn_tf_density["tfoverlap_bin"].divide(syn_tf_density.syn_len)
syn_tf_density = pd.merge(syn_tf_density, df[["syn_id", "core_remodeling_2", "core_mrca_2", "counts"]], how = "left").drop_duplicates()
syn_tf_density.head()
#%%
x, y = "arch", "tf_density"

fig_id = "4A"
re = RE
trim_len = "raw"
frac = "all"
dataset = "all_fantom_syn_non_zero"


plot_figure3(x, y, syn_tf_density.loc[syn_tf_density.tf_density>0], fig_id, re, trim_len, frac, dataset)
#%%
x, y = "arch", "tf_density"

fig_id = "4A"
re = RE
trim_len = "raw"
frac = "all"
dataset = "all_fantom_syn"


plot_figure3(x, y, syn_tf_density, fig_id, re, trim_len, frac, dataset)

#%%

groupby = "arch"
val = "syn_id"
overlap = get_tfbs_overlap_fraction(syn_tf_density, groupby, val)
#%% plot enhancer TF density (including zeros)


fig, ax = plt.subplots(figsize = (6,6))
x = "overallarch"
y = "tf_density"
data = tf_density
outf = make_pdf("all_fantom_enh_x_encode3_tf_density_ENH", RE)

sns.barplot(x, y, data = data, estimator = np.median)

simple_tfden = tf_density.loc[tf_density.overallarch == "simple", "tf_density"]
complex_tfden = tf_density.loc[tf_density.overallarch == "complex", "tf_density"]
print(tf_density.groupby("overallarch")["tf_density"].median())


plt.savefig(outf, bbox_inches = 'tight')

print("simple v. complex", stats.mannwhitneyu(simple_tfden, complex_tfden))
#%% # without zeros


plot_data = tf_density.loc[tf_density.tfoverlap_bin>0]

fig, ax = plt.subplots(figsize = (6,6))
x = "overallarch"
y = "tf_density"
data = plot_data
outf = make_pdf("all_fantom_enh_x_encode3_tf_density_nozero_ENH", RE)

sns.barplot(x, y, data = data, estimator = np.median)

simple_tfden = plot_data.loc[plot_data.overallarch == "simple", "tf_density"]
complex_tfden = plot_data.loc[plot_data.overallarch == "complex", "tf_density"]
print(plot_data.groupby("overallarch")["tf_density"].median())


plt.savefig(outf, bbox_inches = 'tight')

print("simple v. complex", stats.mannwhitneyu(simple_tfden, complex_tfden))


#%% calculate syn TF density


tf_density = df.groupby(["syn_id", "syn_len", "arch"])["tfoverlap_bin"].sum().reset_index()
tf_density["tf_density"] = tf_density["tfoverlap_bin"].divide(tf_density.syn_len)
tf_density.head()

#%% zero overlaps ?

groupby = "arch"
val = "syn_id"
overlap = get_tfbs_overlap_fraction(tf_density, groupby, val)
overlap
#%% plot syntenic block TF density


fig, ax = plt.subplots(figsize = (6,6))
x = "arch"
y = "tf_density"
data = tf_density
outf = make_pdf("all_fantom_enh_x_encode3_tf_density_SYN", RE)

sns.barplot(x, y, data = data, estimator = np.median)
ax.set(xticklabels = ["simple", "core", "derived"])
simple_tfden = tf_density.loc[tf_density.arch == "simple", "tf_density"]
core_tfden = tf_density.loc[tf_density.arch == "complex_core", "tf_density"]
derived_tfden = tf_density.loc[tf_density.arch == "complex_derived", "tf_density"]

print(tf_density.groupby("arch")["tf_density"].median())

print("simple v. core", stats.mannwhitneyu(simple_tfden, core_tfden))
print("core v. derived", stats.mannwhitneyu(derived_tfden, core_tfden))

plt.savefig(outf, bbox_inches = 'tight')



#%%

collection_dict = {}

for tf in df.tf.unique():

    for arch in df.arch.unique():

        comparison_name = tf + "-" + arch
        obs = prep_2x2(tf, arch, df)

        results = quantify_2x2(obs, comparison_name)

        collection_dict[comparison_name] = results
#%%


#%%

results_df = fdr_correction(collection_dict)
#print(len(results_df))
results_df["arch"] = results_df["comparison_name"].apply(lambda x: x.split("-")[1])
results_df["tf"] = results_df["comparison_name"].apply(lambda x: x.split("-")[0])
results_df["log2"]= np.log2(results_df["OR"])

sigresults = results_df.loc[results_df.reject_null == True]

for arch in sigresults.arch.unique():
    print(arch)
    sig_tf = sigresults.loc[sigresults.arch == arch, "tf"].unique()

    plot = sigresults.loc[sigresults.tf.isin(sig_tf)].sort_values(by = ['arch', "log2"])
    x = "tf"
    y = "log2"
    hue = "arch"
    data = plot
    fig, ax = plt.subplots(figsize = (8,24))
    sns.barplot(x=y, y=x, data=data , hue = hue)
    ax.set(xlabel = "OR (log2-scaled)", title = "sig %s" %arch)
    #ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)


#%%


simple_collection_dict = {}

for tf in df.tf.unique():

    arch = "simple"

    comparison_name = tf + "-" + arch
    obs = prep_2x2(tf, arch, df)

    results = quantify_2x2(obs, comparison_name)

    simple_collection_dict[comparison_name] = results
#%%

simple_results_df = fdr_correction(simple_collection_dict)


#%%

simple_results_df.loc[simple_results_df.reject_null == True]
simple_results_df["arch"] = simple_results_df["comparison_name"].apply(lambda x: x.split("-")[1])
simple_results_df["tf"] = simple_results_df["comparison_name"].apply(lambda x: x.split("-")[0])
simple_results_df["log2"]= np.log2(simple_results_df["OR"])

simple_sigresults = simple_results_df.loc[simple_results_df.reject_null == True]

for arch in simple_sigresults.arch.unique():
    print(arch)
    sig_tf = simple_sigresults.loc[simple_sigresults.arch == arch, "tf"].unique()

    plot = simple_sigresults.loc[simple_sigresults.tf.isin(sig_tf)].sort_values(by = ['arch', "log2"])
    x = "tf"
    y = "log2"
    hue = "arch"
    data = plot
    fig, ax = plt.subplots(figsize = (8,24))
    sns.barplot(x=y, y=x, data=data , hue = hue)
    ax.set(xlabel = "OR (log2-scaled)", title = "sig %s" %arch)
    #ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
