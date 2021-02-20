import os
import pandas as pd
from scipy import stats
import statsmodels
import statsmodels.api as sm
import subprocess

FANTOMPATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/all_fantom_enh/ages/"
FANTOMFILE = "syn_breaks_all_fantom_enh_ages.bed"
FANTOM = os.path.join(FANTOMPATH, FANTOMFILE)

ENCODEPATH = "/dors/capra_lab/data/encode/encode3_hg38/TF/"
ENCODEFILE = "trimmed_encRegTfbsClusteredWithCells.liftOver.to.hg19.bed"
ENCODE = os.path.join(ENCODEPATH, ENCODEFILE)

INTERSECTIONPATH = FANTOMPATH
INTERSECTIONFILE = "All_FANTOM_x_ENCODE_hg19.bed"
INTERSECTION = os.path.join(INTERSECTIONPATH, INTERSECTIONFILE)

RE = "/dors/capra_lab/projects/enhancer_ages/fantom/results/encode3/"

#%% Functions

def bed_intersect(fantom, encode, intersection):

    cmd = "bedtools intersect -a %s -b %s -wao > %s" % (fantom, encode, intersection)

    subprocess.call(cmd, shell = True)

    print(cmd)


def format_df(df):

    # remove chip-peaks with overlap less than 6bp

    df = df.loc[(df.overlap >5) | (df.overlap ==0)].copy() # removes 34802 TF overlaps


    # add architecture label
    df["arch"] = "simple"
    df.loc[(df.core_remodeling ==1) & (df.core ==1), "arch"] = "complex_core"
    df.loc[(df.core_remodeling ==1) & (df.core ==0), "arch"] = "complex_derived"

    df["overallarch"] = "simple"
    df.loc[df.core_remodeling ==1, "overallarch"] = "complex"
    # add syn identifier
    df["syn_id"] = df.chr_syn + ":" + df.start_syn.map(str) + "-" + df.end_syn.map(str)

    #calculate enhancer and syntenic block length.

    df["enh_len"] = df.end - df.start
    df["syn_len"] = df.end_syn - df.start_syn

    # binary for TF overlap
    df["tfoverlap_bin"] = 1
    df.loc[df.tf == ".", "tfoverlap_bin"] = 0

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


#%% Bed command

bed_intersect(FANTOM, ENCODE, INTERSECTION)

#%% dataframe

cols = ["chr_syn", "start_syn", "end_syn",
"enh_id","chr", "start", "end",
"seg_index", "core_remodeling", "core",
"mrca",
"chr_tf", "start_tf", "end_tf",
"tf_id", "peak_len", "reads",
"tf",  "overlap"
]

df_ = pd.read_csv(INTERSECTION,
sep = '\t',
header = None)

df_.columns = cols # add column names

df_.info()

df = format_df(df_) # format the dataframe

df.info()

df.head()
#%%

#%%
df.loc[df.enh_id == "chr9:99132173-99132333"][["mrca", "core", "tf", "syn_len", "tfoverlap_bin"]]
#%% get some basic info about Fantom enhancer overlap


enh_df = df.groupby(["enh_id", "core_remodeling", "overallarch"])[["mrca", "seg_index"]].max().reset_index()

totalenh_n = len(enh_df) #30279 enhancers total
simpleenh_n = len(enh_df.loc[enh_df.overallarch == "simple"]) #14098 simple enhancers
complexenh_n = len(enh_df.loc[enh_df.overallarch != "simple"]) # 8744 complex enhancers


#%% calculate enhancer TF density


tf_density = df.groupby(["enh_id", "enh_len", "overallarch"])["tfoverlap_bin"].sum().reset_index().drop_duplicates()

tf_density["tf_density"] = tf_density["tfoverlap_bin"].divide(tf_density.enh_len)



#%% how many enhancers do not overlap  TFs?


tf_density.loc[tf_density.tfoverlap_bin == 0].groupby("overallarch")["enh_id"].count()

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
tf_density.loc[tf_density.tfoverlap_bin == 0].groupby("arch")["syn_id"].count()

core_n = tf_density.loc[tf_density.arch == "complex_core"]["syn_id"].count()
derived_n = tf_density.loc[tf_density.arch == "complex_derived"]["syn_id"].count()

print(core_n, derived_n)

4031/11793  # cores with no overlaps
6820/15749 # derived with no overlaps

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
