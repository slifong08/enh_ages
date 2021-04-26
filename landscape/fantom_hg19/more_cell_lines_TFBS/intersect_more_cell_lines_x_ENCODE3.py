import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

from scipy import stats
import seaborn as sns
import statsmodels
import statsmodels.api as sm
import subprocess

FANTOMBASE = "/dors/capra_lab/projects/enhancer_ages/fantom/data/ENCODE3/"

ENCODEPATH = "/dors/capra_lab/data/encode/encode3_hg38/TF/liftOver_hg19/"

RE = "/dors/capra_lab/projects/enhancer_ages/landscape/results/fantom_x_encode3_tfbs/"

colors = [ "amber", "dusty purple", "windows blue"]
PAL = sns.xkcd_palette(colors)
sns.palplot(PAL)

colors = [ "windows blue"]
DERPAL = sns.xkcd_palette(colors)
sns.palplot(DERPAL)


#%% Functions


def get_cell_lines():
    sample_dict = {
    "all_fantom_enh": "all_fantom_enh",
    "A549":"A549_FANTOM5_tpm_hg19",
    "GM12878_CL":"CL_0000945_lymphocyte_of_B_lineage_expressed_enhancers",
    "GM12878": "GM12878_FANTOM5_tpm_hg19",
    "HepG2": "HEPG2_FANTOM5_tpm_hg19",
    "HepG2_CL": "CL_0000182_hepatocyte_expressed_enhancers",
    "K562":"CL_0000094_granulocyte_expressed_enhancers",
    "liver":"UBERON_0002107_liver_expressed_enhancers",
    #"H1-hESC":"H1-esc_FANTOM5_tpm_hg19",
    #"MCF7":"CL_0002327_mammary_epithelial_cell_expressed_enhancers",
    #"PC-3": "PC-3_FANTOM5_tpm_hg19",
    }

    return sample_dict


def get_paths(cell_line, file_tag, fantombase, encodepath):

    FANTOMPATH = os.path.join(fantombase, file_tag, "ages")
    FANTOMFILE = "syn_breaks_%s_ages.bed" % file_tag
    FANTOM = os.path.join(FANTOMPATH, FANTOMFILE)

    if "CL" in cell_line:
        ENCODEFILE = "cells/" + cell_line.split("_CL")[0] + ".bed.gz"
    elif cell_line == "all_fantom_enh":
        ENCODEFILE = "trimmed_encRegTfbsClusteredWithCells.liftOver.to.hg19.bed"
    else:
        ENCODEFILE = "cells/" + cell_line + ".bed.gz"

    ENCODE = os.path.join(encodepath, ENCODEFILE)


    INTERSECTIONPATH = os.path.join(fantombase, file_tag)
    INTERSECTIONFILE = "%s_x_ENCODE.bed" % file_tag
    INTERSECTION = os.path.join(INTERSECTIONPATH, INTERSECTIONFILE)

    return FANTOM, ENCODE, INTERSECTION


def bed_intersect(fantom, encode, intersection):

    if os.path.exists(intersection) == False:
        cmd = "bedtools intersect -a %s -b %s -wao > %s" % (fantom, encode, intersection)

        subprocess.call(cmd, shell = True)

        print(cmd)
    else:
        print("you have previously done the enh x encode intersection")


def format_df(intersection_file):

    cols = ["chr_syn", "start_syn", "end_syn",
    "enh_id","chr", "start", "end",
    "seg_index", "core_remodeling", "core",
    "mrca",
    "chr_tf", "start_tf", "end_tf",
    "tf_id", "peak_len",
    "tf", "cell_line" , "overlap"
    ]

    df = pd.read_csv(intersection_file,
    sep = '\t',
    header = None).drop_duplicates()

    df.columns = cols # add column names

    df["tf"] = df["tf_id"].apply(lambda x: x.split("_")[0])

    # add architecture label - core, derived, simple
    df["arch"] = "simple"
    df.loc[(df.core_remodeling ==1) & (df.core ==1), "arch"] = "complex_core"
    df.loc[(df.core_remodeling ==1) & (df.core ==0), "arch"] = "complex_derived"

    # add architecture label - complex, simple
    df["overallarch"] = "simple"
    df.loc[df.core_remodeling ==1, "overallarch"] = "complex"

    # add syn identifier
    df["syn_id"] = df.chr_syn + ":" + df.start_syn.map(str) + "-" + df.end_syn.map(str)

    #calculate enhancer and syntenic block length
    df["enh_len"] = df.end - df.start
    df["syn_len"] = df.end_syn - df.start_syn
    df.loc[df.syn_len <6, "syn_len"] = 0


    # binary for TF overlap, any TF that overlaps less than 6bp is not counted.
    df["tfoverlap_bin"] = 1
    df.loc[df.tf == ".", "tfoverlap_bin"] = 0
    df.loc[df.overlap <6 , "tfoverlap_bin"] = 0

    return df


def count_enhancers(df, arch):

    if arch == "enh":
        enh_df = df.groupby(["enh_id", "core_remodeling", "overallarch"])[["mrca", "seg_index"]].max().reset_index()

        totalenh_n = len(enh_df) #30279 enhancers total
        simpleenh_n = len(enh_df.loc[enh_df.overallarch == "simple"]) #14098 simple enhancers
        complexenh_n = len(enh_df.loc[enh_df.overallarch != "simple"]) # 8744 complex enhancers


        return totalenh_n, simpleenh_n, complexenh_n

    elif arch == "syn":

        total = len(df)
        core_n = df.loc[df.arch == "complex_core"]["syn_id"].count()
        derived_n = df.loc[df.arch == "complex_derived"]["syn_id"].count()
        simple_n = df.loc[df.arch == "simple"]["syn_id"].count()

        return total, core_n, derived_n, simple_n


def mwu(tf_density, arch):

    # calculate means
    median = (tf_density.groupby("arch")["tf_density"].median())
    print("\narchitecture medians", median)

    if arch == "enh":

        #stratify dataframe by simple and complex arch
        simple_tfden = tf_density.loc[tf_density.arch == "simple", "tf_density"]
        complex_tfden = tf_density.loc[tf_density.arch == "complex", "tf_density"]

        # calculate MWU
        test_arch, p_arch = stats.mannwhitneyu(simple_tfden, complex_tfden)
        print("\n", "simple v. complex enh MWU stat =", round(test_arch,3), "p =", p_arch )

        return test_arch, p_arch, median

    elif arch == "syn":

        simple_tfden = tf_density.loc[tf_density.arch == "simple", "tf_density"]
        core_tfden = tf_density.loc[tf_density.arch == "complex_core", "tf_density"]
        derived_tfden = tf_density.loc[tf_density.arch == "complex_derived", "tf_density"]

        testcore, pcore = stats.mannwhitneyu(simple_tfden, core_tfden)
        print("\n", "simple v. complexcore MWU stat =", round(testcore,3), "p =", pcore)

        test_der, p_der = stats.mannwhitneyu(derived_tfden, core_tfden)
        print("\n", "core v. derived MWU stat =", round(test_der,3), "p =", p_der)

        return testcore, pcore, test_der, p_der, median


def calculate_tf_density(arch, df):

    density_cols = ["id", "len", "arch", "tfoverlap_bin", "tf_density",]


    if arch == "enh":

        tf_density = df.groupby(["enh_id", "enh_len", "overallarch"])["tfoverlap_bin"].sum().reset_index().drop_duplicates()
        tf_density["tf_density"] = tf_density["tfoverlap_bin"].divide(tf_density.enh_len)


    elif arch == "syn":

        tf_density = df.groupby(["syn_id", "syn_len", "arch"])["tfoverlap_bin"].sum().reset_index()
        tf_density["tf_density"] = tf_density["tfoverlap_bin"].divide(tf_density.syn_len)

    # rename columns
    tf_density.columns = density_cols

    # how many enhancers do not overlap  TFs?
    zero_overlap = tf_density.loc[tf_density.tfoverlap_bin == 0].groupby("arch")["id"].count().reset_index()

    return tf_density, zero_overlap


def calculate_zero_syn_freq(zero_syn, df, cell_line, RE):

    zero_syn.columns = ['arch', "zero_counts"]

    arch_df = df[["arch", "syn_id"]].drop_duplicates()

    total_arch_counts = arch_df.groupby(["arch"])["syn_id"].count().reset_index()
    total_arch_counts.columns = ['arch', "total_counts"]

    zero_syn = pd.merge(zero_syn, total_arch_counts, on = "arch")
    zero_syn["freq_zero"] = zero_syn.zero_counts.divide(zero_syn.total_counts)
    zero_syn["freq_nonzero"] = 1-zero_syn["freq_zero"]

    print(zero_syn)
    zero_syn.to_csv('%snonzero_%s.csv' % (RE, cell_line), index = False)


def plot_bar_tf_density(x, y, data, outf, order, p, med):

    # plot enhancer TF density
    fig, ax = plt.subplots(figsize = (6,6))
    sns.set("poster")

    sns.barplot(x, y, data = data, estimator = np.median,
    order = order, palette = PAL, n_boot = 10000)

    ax.set(
    xlabel= "p = %s\n%s" % (p,med),
    title = outf.split("/")[-1],
    ylabel = "TFBS density\nmedian"
    )

    plt.savefig(outf, bbox_inches = 'tight')


def make_pdf(file_name, RE):

    OUTFILE = file_name + ".pdf"
    OUTF = os.path.join(RE, OUTFILE)

    return OUTF


def prep_2x2(tf, arch1, arch2, df):

    # only evaluate enhancers that overlap TF.
    # This excludes enhancers with zero overlaps from background set.
    df = df.loc[df.tfoverlap_bin >0]

    # split dataframe by two architectures to compare
    dfarch = df.loc[df.arch == arch1]

    if arch2 == "bkgd":
        arch2 = "all_enh_bkgd"
        dfbkgd = df.loc[df.arch != arch1]

    else:
        dfbkgd = df.loc[df.arch == arch2]

    comparison_name = tf + "-" + arch1 + "_v_" + arch2

    # count how many TF overlaps are in each arch.
    TF_in_arch = len(dfarch.loc[dfarch.tf == tf])
    TF_bkgd = len(dfbkgd.loc[dfbkgd.tf == tf])
    not_TF_in_arch = len(dfarch.loc[dfarch.tf != tf])
    not_TF_in_bkgd = len(dfbkgd.loc[dfbkgd.tf != tf])

    a, b, c, d = TF_in_arch, not_TF_in_arch, TF_bkgd, not_TF_in_bkgd

    obs = [[a,b], [c,d]]

    return obs, comparison_name


def quantify_2x2(obs, comparison_name, min_instances):

    if obs[0][0] > min_instances or obs[1][0]>min_instances:

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

        return newdf


def fdr_correction(collection_dict, alpha):

    df = pd.concat(collection_dict.values())

    pvals = df["P"]

    df["reject_null"], df["FDR_P"] = statsmodels.stats.multitest.fdrcorrection(pvals, alpha=alpha)

    # other dataframe formatting
    df["arch"] = df["comparison_name"].apply(lambda x: x.split("-")[1])
    df["tf"] = df["comparison_name"].apply(lambda x: x.split("-")[0])
    df["log2"]= np.log2(df["OR"])

    return df


def plot_bar_tf_enrichment(df, cell_line, outf, alpha):

    fig, ax = plt.subplots(figsize = (6,12))
    sns.set("poster")

    x = "tf"
    y = "log2"
    hue = "arch"
    data = df.sort_values(by = y)

    sns.barplot(x=y, y=x, data=data , hue = hue, palette = DERPAL)

    ax.legend(bbox_to_anchor = (1,1))
    ax.set(xlabel = "OR log2-scale\n FDR<%s" % str(alpha), title = cell_line)

    plt.savefig(outf, bbox_inches = "tight")


def run_2x2(arch1, arch2, df, min_instances, alpha):

    collection_dict = {}

    for tf in df.tf.unique():

        if tf != ".":

            obs, comparison_name = prep_2x2(tf, arch1, arch2, df)

            results = quantify_2x2(obs, comparison_name, min_instances)

            collection_dict[comparison_name] = results

    #FDR correction

    results_df = fdr_correction(collection_dict, alpha)

    df = results_df.loc[results_df.reject_null == True]

    if len(df)>0:
        outf = make_pdf("%s_enh_x_encode3_sig_tf_arch_enrichment_%s_v_%s_FDR_%s" % (cell_line, arch1, arch2, alpha), RE)

        plot_bar_tf_enrichment(df, cell_line, outf, alpha)

    else:
        print("\nno sig results for comparison", arch1, "v.", arch2)

    return results_df


def run_analysis(cell_line, val, fantombase, encodepath, min_instances, alpha):

    print(cell_line, val)
    fantom, encode, intersection = get_paths(cell_line, val, fantombase, encodepath)

    #Bed command
    bed_intersect(fantom, encode, intersection)

    #dataframe
    df = format_df(intersection)

    #get some basic info about Fantom enhancer overlap
    arch = "enh"
    totaln, simplen, complexn = count_enhancers(df, arch)

    # calculate enhancer TF density
    tf_density_enh, zero_enh = calculate_tf_density(arch, df)


    # plot all enhancer-level data
    x, y = "arch", "tf_density"
    order = ["simple", "complex"]

    data = tf_density_enh
    outf = make_pdf("%s_enh_x_encode3_tf_density_%s"  % (cell_line, arch), RE)

    test_arch, p_arch, median = mwu(tf_density_enh, arch)
    plot_bar_tf_density(x, y, data, outf, order, p_arch, median)


    print("\nNon-zero TFBS densities only")


    # plot all enhancer-level data without zeros
    non_zero_tf_density = tf_density_enh.loc[tf_density_enh.tfoverlap_bin>0]

    data = non_zero_tf_density
    outf = make_pdf("%s_enh_x_encode3_tf_density_%s_non_zero_tf_density" % (cell_line, arch), RE)
    test_arch_, p_arch_, median = mwu(non_zero_tf_density, arch)
    plot_bar_tf_density(x, y, data, outf, order, p_arch_, median)



    # calculate syn-level TF density
    arch = "syn"
    totaln, coren, derivedn, simplen = count_enhancers(df, arch)
    tf_density_syn, zero_syn = calculate_tf_density(arch, df)

    # calculate frequency of derived sequences that do not overlap TFBS

    calculate_zero_syn_freq(zero_syn, df, cell_line, RE)

    print("\nSyn TFBS densities")

    # plot syn block TF density
    order = ["simple", "complex_core", "complex_derived"]
    data = tf_density_syn
    outf = make_pdf("%s_enh_x_encode3_tf_density_%s" % (cell_line, arch), RE)

    testcore, pcore, test_der, p_der, median = mwu(tf_density_syn, arch)
    new_p = "simple_v_core p = %s,  core_v_der= %s" %(pcore, p_der)
    plot_bar_tf_density(x, y, data, outf, order, new_p, median)


    print("\nNon-zero syn TFBS densities only")
    non_zero_syn_tf_density = tf_density_syn.loc[tf_density_syn.tfoverlap_bin>0]


    # plot non-zero syn block TF density
    data = non_zero_syn_tf_density
    outf = make_pdf("%s_syn_x_encode3_tf_density_%s_non_zero_tf_density" % (cell_line, arch), RE)

    testcore, pcore, test_der, p_der, median = mwu(non_zero_syn_tf_density, arch)
    new_p = "simple_v_core p = %s,  core_v_der= %s" %(pcore, p_der)
    plot_bar_tf_density(x, y, data, outf, order, new_p, median)


    ### TF ENRICHMENT IN ARCHITECTURE ###
    # DER V. CORE
    # DER V. BKGD

    # calculate TF enrichment in architecture/syn blocks
    arch1, arch2 = "complex_derived", "complex_core"
    der_v_core = run_2x2(arch1, arch2, df, MIN_INSTANCES, ALPHA)

    arch1, arch2 = "complex_derived", "bkgd"
    der_v_bkgd = run_2x2(arch1, arch2, df, MIN_INSTANCES, ALPHA)

    return der_v_core, der_v_bkgd, tf_density_enh, tf_density_syn, df


#%%
sample_dict = get_cell_lines()

results_dict = {}
der_v_core_dict, der_v_bkgd_dict = {}, {} # collect all the dataframes for tf enrichment
tf_den_enh = {}
tf_den_syn = {}

#%%

ALPHA = 0.05
MIN_INSTANCES = 500


cell_line = "all_fantom_enh"
val = sample_dict[cell_line]


der_v_core, der_v_bkgd, tf_density_enh, tf_density_syn, df = run_analysis(cell_line, val, FANTOMBASE, ENCODEPATH, MIN_INSTANCES, ALPHA)

results_dict[cell_line] = df
tf_den_enh[cell_line] = tf_density_enh
tf_den_syn[cell_line] = tf_density_syn
der_v_bkgd_dict[cell_line] = der_v_bkgd
der_v_core_dict[cell_line] = der_v_core

der_v_core.head()

#%%
outf = f"{RE}derived_v_core_OR_all_fantom_x_encode3_tfbs.tsv"
der_v_core.to_csv(outf, sep = '\t', index = False)

outf = f"{RE}derived_v_bkgd_OR_all_fantom_x_encode3_tfbs.tsv"
der_v_bkgd.to_csv(outf, sep = '\t', index = False)

#%% run enrichment per cell_line

ALPHA = 0.1
MIN_INSTANCES = 50

cell_lines = ["HepG2", "HepG2_CL", "K562", "GM12878", "GM12878_CL",  "A549", "liver"]
for cell_line in cell_lines:
    val = sample_dict[cell_line]


    der_v_core, der_v_bkgd, tf_density_enh, tf_density_syn, df = run_analysis(cell_line, val, FANTOMBASE, ENCODEPATH, MIN_INSTANCES, ALPHA)

    results_dict[cell_line] = df
    tf_den_enh[cell_line] = tf_density_enh
    tf_den_syn[cell_line] = tf_density_syn
    der_v_bkgd_dict[cell_line] = der_v_bkgd
    der_v_core_dict[cell_line] = der_v_core


#%%

df = tf_den_syn["all_fantom_enh"]
fulldf = results_dict["all_fantom_enh"]
syn_ages = fulldf[["syn_id", "mrca", "enh_id"]].drop_duplicates()
syn_ages.columns= ["id", "mrca", "enh_id"]
df = pd.merge(df, syn_ages)

core_age = fulldf.groupby(["enh_id"])["mrca"].max().reset_index()
core_age.columns= ["enh_id", 'core_mrca']
df = pd.merge(df, core_age)

df.head()
#%%

# age and taxon file
syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2"]] # whittle down the df
df["core_mrca"] = df["core_mrca"].round(3) # round the ages

df = pd.merge(df, syn_gen_bkgd, how = "left", left_on = 'core_mrca', right_on = "mrca")
#%%
df.groupby(["core_mrca", "arch"])["tf_density"].median()
#%%
x = "mrca_2"
y = "tf_density"
hue = "arch"
data = df
xlabs = ["homo", "prim", "euar", "bore", "euth", "ther", "mam", "amni", "tetr", "vert"]

fig, ax = plt.subplots(figsize = (10,10))
sns.barplot(x = x, y = y, data = data, hue = hue, palette = PAL, estimator = np.median)#order = ["simple", "complex_core", "complex_derived"])

ax.set(ylabel = "tfbs density\nmedian", title = "all_fantom_enh", xlabel = "core_age")
ax.set_xticklabels(xlabs, rotation = 90)
ax.legend(bbox_to_anchor = (1,1))

plt.savefig("%sall_fantom_tfbs_mrca.pdf"%RE, bbox_inches="tight")


#%%
x = "mrca_2"
y = "tf_density"
hue = "arch"
data = df.loc[df[y]>0]

fig, ax = plt.subplots(figsize = (10,10))
sns.barplot(x = x, y = y, data = data, hue = hue, palette = PAL, estimator = np.median)#order = ["simple", "complex_core", "complex_derived"])
xlabs = ["homo", "prim", "euar", "bore", "euth", "ther", "mam", "amni", "tetr", "vert"]
ax.set(ylabel = "tfbs density\nmedian", title = "all_fantom_enh", xlabel = "core_age")
ax.set_xticklabels(xlabs, rotation = 90)
ax.legend(bbox_to_anchor = (1,1))

plt.savefig("%sall_fantom_syn_non_zero_tfbs_mrca.pdf"%RE, bbox_inches="tight")
#%%

df["core"] = "core"
df.loc[df.arch == "complex_derived", "core"] = "der"
df.loc[df.len <6, "len"] = 0


complex_only = df.loc[df.arch.str.contains("complex") & (df.len > 0)]


regions = complex_only.groupby(["enh_id", "core"])[["tfoverlap_bin", "len"]].sum().reset_index() # summarize complex enhancers by regions
regions["tf_density"] = regions.tfoverlap_bin.divide(regions.len) # calculate the total regional syn density

regions.head()
#%%

coreder_den = regions.pivot(index = "enh_id", columns = "core", values = "tf_density").reset_index()
coreder_len = regions.pivot(index = "enh_id", columns = "core", values = "len").reset_index()

coreder_len["sum_len"] = coreder_len.der + coreder_len.core # get the full length of the enhancer
coreder_len["ratio_der"] = coreder_len.der.divide(coreder_len["sum_len"]) # get the ratio of derived length per enhancer
coreder_len["pct"] = coreder_len["sum_len"].rank(pct = True) # rank lengths

coreder_len.head()

coreder_len.describe()
#%%
'''

core	core	der	sum_len	ratio_der	pct
count	10758.000000	10725.000000	10528.000000	10528.000000	10528.000000
mean	199.353783	177.452867	373.782390	0.461647	0.500047
std	145.930430	158.158972	187.829319	0.281999	0.288688
min	6.000000	6.000000	15.000000	0.007609	0.000095
25%	89.000000	62.000000	257.000000	0.209934	0.250237
50%	175.000000	141.000000	350.000000	0.442009	0.500285
75%	278.000000	249.000000	436.000000	0.703121	0.749430
max	1669.000000	2791.000000	2860.000000	0.993852	1.000000
'''
#%% assess enhancer length  v. derived ratio inner quartiles
iqrlen = coreder_len.loc[(coreder_len.pct<= 0.75) & (coreder_len.pct>=0.25)]
ids = iqrlen.enh_id
iqrlen.shape
x = "sum_len"
y = "ratio_der"

outf = f"{RE}all_fantom_enh_lenIRQ_x_ratio_der.pdf"
data = iqrlen
sns.jointplot(x=x, y = y, data = data, kind = "hex")
iqrlen.shape
iqrlen.describe()

plt.savefig(outf, bbox_inches = "tight")
#%%
iqrden = coreder_den.loc[coreder_den.enh_id.isin(iqrlen.enh_id)]
iqrden.head()

#%%
iqrden['total'] = iqrden.core +  iqrden.der
iqrden = iqrden.loc[(iqrden['core'] >0) & (iqrden['der'] >0)]
iqrden.shape

iqrden["ratio_der_core_den"] = iqrden.der/iqrden.core
iqrden.describe()


#%%
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(
         iqrden.iloc[:, 1], iqrden.iloc[:, 2], random_state=0)
X_train = np.array(X_train).reshape(-1,1)
y_train = np.array(y_train).reshape(-1,1)


reg = LinearRegression().fit(X_train, y_train)
reg.score(X_train, y_train)
#%%
len(X_train)
#%%
x = X_train[:, 0] # derived
y = y_train[:, 0] # core

sns.jointplot(x=x, y = y,)

outf = f"{RE}all_fantom_enh_lenIRQ_core_v_der_tfbs_den.pdf"

plt.savefig(outf, bbox_inches = "tight")
#%%
sns.histplot(iqrden.loc[iqrden.ratio_der_core_den<1.76, "ratio_der_core_den"])
#%%
