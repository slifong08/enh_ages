import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy import stats
import seaborn as sns
import statsmodels
import statsmodels.api as sm
import subprocess

ENHBASE = "/dors/capra_lab/projects/enhancer_ages/encode/data/"
#ENHBASE = "/dors/capra_lab/projects/enhancer_ages/encode/hepg2/data/"

ENCODEPATH = "/dors/capra_lab/data/encode/encode3_hg38/TF/"

RE = "/dors/capra_lab/projects/enhancer_ages/landscape/results/cCRE_x_tfbs_encode3/HepG2/"
RE_DATA = RE + "data/"

if os.path.exists(RE_DATA) == False:
    os.mkdir(RE_DATA)

colors = ["amber", "dusty purple", "windows blue"]
PAL = sns.xkcd_palette(colors)
sns.palplot(PAL)

colors = ["windows blue"]
DERPAL = sns.xkcd_palette(colors)
sns.palplot(DERPAL)


#%% Functions

def get_cell_lines():
    sample_dict = {
        "HepG2": "ELS_combined_HepG2",
    }
    return sample_dict


def get_paths(cell_line, file_tag, fantombase, encodepath):

    FANTOMPATH = os.path.join(fantombase, file_tag, "ages")
    FANTOM = os.path.join(FANTOMPATH, f"{file_tag}syn_breaks_%s_ages.bed")

    if "CL" in cell_line:
        ENCODEFILE = "cells/" + cell_line.split("_CL")[0] + ".bed.gz"
    elif cell_line == "all_fantom_enh":
        ENCODEFILE = "trimmed_encRegTfbsClusteredWithCells.bed"
    else:
        ENCODEFILE = "cells/" + cell_line + ".bed.gz"

    ENCODE = os.path.join(encodepath, ENCODEFILE)

    INTERSECTIONPATH = os.path.join(fantombase, file_tag)
    INTERSECTION = os.path.join(INTERSECTIONPATH, f"{file_tag}_x_ENCODE.bed")

    return FANTOM, ENCODE, INTERSECTION


def bed_intersect(fantom, encode, intersection):

    if os.path.exists(intersection) is False:
        cmd = f"bedtools intersect -a {fantom} -b {encode} -wao > {intersection}"

        subprocess.call(cmd, shell=True)

        print(cmd)
    else:
        print("previously done enh x encode intersection")


def get_core_age(df):
    core = df.groupby("enh_id")["mrca_2"].max().reset_index()
    core.columns = ["enh_id", 'core_mrca_2']
    df = pd.merge(df, core, how="left")

    return df


def reEval_PrimComplex(enh):

    # get all the complex enhancers w/ primate core ages
    prComEnhID = enh.loc[(enh.core == 1)
                         & (enh.core_remodeling == 1)
                         & (enh.taxon2.str.contains("Primate"))]["enh_id"].unique()

    # get all the complex enhancer ids where there is a real human derived sequence
    pr_complex = enh.loc[(enh.enh_id.isin(prComEnhID))
                         & (enh.core_remodeling == 1)
                         & (enh.core == 0)
                         & (enh.mrca == 0),
                         ]["enh_id"]

    # i'm going to reassign any primate complex enhancer
    # where derived regions are from other primates
    # get the set of primate complex enhancers w/ primate derived sequences
    # and rename them as simple enhancers
    pr_simple = set(prComEnhID) - set(pr_complex)

    # reassign core and core remodeling columns
    enh.loc[enh.enh_id.isin(pr_simple), "core"] = 1
    enh.loc[enh.enh_id.isin(pr_simple), "core_remodeling"] = 0
    return enh


def reassign_small_taxons(df):
    # lump small samples with larger, older ancestors
    df.loc[df.taxon2 == "Sarcopterygian", "taxon2"] = "Vertebrata"
    df.loc[df.taxon2 == "Tetrapoda", "taxon2"] = "Vertebrata"
    df.loc[df.taxon2 == "Euarchontoglires", "taxon2"] = "Boreoeutheria"
    df.loc[df.taxon2 == "Primates", "taxon2"] = "Boreoeutheria"
    return df


def merge_syn_annot(df):
    SYN_GROUP = "/dors/capra_lab/projects/enhancer_ages/hg38_syn_taxon.bed"
    syn = pd.read_csv(SYN_GROUP, sep='\t')

    # round all values

    syn[["mrca", "mrca_2"]] = syn[["mrca", "mrca_2"]].round(3)
    df["mrca"] = df["mrca"].round(3)

    df = pd.merge(df, syn, how="left")

    return df


def format_df(intersection_file):

    # do the dance - add mrca_2 column, get mrca_2 core age, drop mrca_2 column, then add it back, but this time to reflect the core_age and core taxon, instead of the syntenic age.

    cols = ["chr_syn", "start_syn", "end_syn",
            "enh_id", "chr", "start", "end",
            "seg_index", "core_remodeling", "core",
            "mrca",
            "chr_tf", "start_tf", "end_tf",
            "tf_id", "peak_len",
            "tf", "cell_line", "overlap"
            ]

    df = pd.read_csv(intersection_file,
                     sep='\t',
                     header=None).drop_duplicates()

    df.columns = cols  # add column names

    df = merge_syn_annot(df) # merge taxon information

    df["tf"] = df["tf_id"].apply(lambda x: x.split("_")[0])

    # add architecture label - core, derived, simple
    df["arch"] ="complex_core"
    df.loc[(df.core_remodeling == 0), "arch"] =  "simple"
    df.loc[df.core == 0, "arch"] = "complex_derived"

    df = reEval_PrimComplex(df)

    # add architecture label - complex, simple
    df["overallarch"] = "simple"
    df.loc[df.core_remodeling == 1, "overallarch"] = "complex"

    # add syn identifier
    df["syn_id"] = df.chr_syn + ":" + \
        df.start_syn.map(str) + "-" + df.end_syn.map(str)

    #calculate enhancer and syntenic block length
    df["enh_len"] = df.end - df.start
    df["syn_len"] = df.end_syn - df.start_syn
    df.loc[df.syn_len < 6, "syn_len"] = 0

    # binary for TF overlap, any TF that overlaps less than 6bp is not counted.
    df["tfoverlap_bin"] = 1
    df.loc[df.tf == ".", "tfoverlap_bin"] = 0
    df.loc[df.overlap < 6, "tfoverlap_bin"] = 0

    return df


def count_enhancers(df, arch, stats_f):

    if arch == "enh":
        enh_df = df.groupby(["enh_id", "core_remodeling", "overallarch"])[
                            ["mrca", "seg_index"]].max().reset_index()

        totalenh_n = enh_df.shape[0]  # 30279 enhancers total
        # 14098 simple enhancers
        simpleenh_n = enh_df.loc[enh_df.overallarch == "simple"].shape[0]
        # 8744 complex enhancers
        complexenh_n = enh_df.loc[enh_df.overallarch != "simple"].shape[0]
        info = f"total enh N - {totalenh_n}\n simple enh N - {simpleenh_n}\n\
         complex enhancer N - {complexenh_n}"


    elif arch == "syn":

        total = df.shape[0]
        core_n = df.loc[df.arch == "complex_core"]["syn_id"].count()
        derived_n = df.loc[df.arch == "complex_derived"]["syn_id"].count()
        simple_n = df.loc[df.arch == "simple"]["syn_id"].count()

        pt1 = f"total syn N - {total}"
        pt2 = f"\nsimple syn N - {simple_n}"
        pt3 = f"\ncomplex core N - {core_n}"
        pt4 = f"\ncomplex derived N - {derived_n}"
        info = pt1 + pt2 + pt3 +pt4

    with open(stats_f, "a") as write_f:
        write_f.write(info)


def mwu(tf_density, arch, kw, stats_f):

    # calculate means
    median = (tf_density.groupby("arch")["tf_density"].median())
    print("\narchitecture medians", median)

    if arch == "enh":

        #stratify dataframe by simple and complex arch
        simple_tfden = tf_density.loc[tf_density.arch
                                      == "simple", "tf_density"]
        complex_tfden = tf_density.loc[tf_density.arch
                                       == "complex", "tf_density"]

        # calculate MWU
        test_arch, p_arch = stats.mannwhitneyu(simple_tfden, complex_tfden)
        p1 = f"\n###{kw}###"
        p2 = "\nsimple v. complex enh TFBS MWU"
        p3 = f"\nstat = {round(test_arch, 3)}, \np = {p_arch}"
        info  = p1 + p2 + p3



        listP = [p_arch]
    elif arch == "syn":

        simple_tfden = tf_density.loc[tf_density.arch\
                                      == "simple", "tf_density"]
        core_tfden = tf_density.loc[tf_density.arch\
                                    == "complex_core", "tf_density"]
        derived_tfden = tf_density.loc[tf_density.arch\
                                       == "complex_derived", "tf_density"]

        testcore, pcore = stats.mannwhitneyu(simple_tfden, core_tfden)

        test_der, p_der = stats.mannwhitneyu(derived_tfden, core_tfden)

        p1 = f"\n###{kw}###"
        p2 = "\nsimple v. complexcore TFBS MWU"
        p3 = f"\nstat = {round(testcore, 3)}\np = {pcore}"
        p4 = "\ncore v. derived MWU"
        p5 = f"\nstat = {round(test_der, 3)},\np = {p_der}"
        info  = p1 + p2 + p3 + p4 + p5

        listP = [pcore, p_der]

    with open(stats_f, "a") as write_f:
        write_f.write(info)

    return listP, median

def calculate_tf_density(arch, df):

    density_cols = ["id", "len", "arch", "tfoverlap_bin", "tf_density", ]

    if arch == "enh":
        cols = ["enh_id", "enh_len", "overallarch"]
        tf_density = df.groupby(cols)["tfoverlap_bin"].sum().reset_index().drop_duplicates()

        tf_density["tf_density"] = tf_density["tfoverlap_bin"].divide(tf_density.enh_len)

    elif arch == "syn":
        cols = ["syn_id", "syn_len", "arch"]
        tf_density = df.groupby(cols)["tfoverlap_bin"].sum().reset_index()
        tf_density["tf_density"] = tf_density["tfoverlap_bin"].divide(tf_density.syn_len)

    # rename columns
    tf_density.columns = density_cols

    # how many enhancers do not overlap  TFs?
    zero_overlap = tf_density.loc[
    tf_density.tfoverlap_bin == 0].groupby("arch")[
    "id"].count().reset_index()

    return tf_density, zero_overlap


def calculate_zero_syn_freq(zero_syn, df, stats_f):

    zero_syn.columns = ['arch', "zero_counts"]

    arch_df = df[["arch", "syn_id"]].drop_duplicates()

    total_arch_counts = arch_df.groupby(
        ["arch"])["syn_id"].count().reset_index()  # counts per arch
    total_arch_counts.columns = ['arch', "total_counts"]

    zero_syn = pd.merge(zero_syn, total_arch_counts, on="arch")
    zero_syn["freq_zero"] = zero_syn.zero_counts.divide(zero_syn.total_counts)
    zero_syn["freq_nonzero"] = 1-zero_syn["freq_zero"]

    print(zero_syn)

    with open(stats_f, "a") as write_f:
        write_f.write(" ".join(list(zero_syn.columns)))

        for index, row in zero_syn.iterrows():
            write_f.write(row.to_string())


def plot_bar_tf_density(x, y, data, outf, order, p, med):

    # plot enhancer TF density
    fig, ax = plt.subplots(figsize=(6, 6))
    sns.set("poster")

    sns.barplot(x, y, data=data, estimator=np.median,
                order=order, palette=PAL, n_boot=10000)

    ax.set(
        xlabel="p = %s\n%s" % (p, med),
        title=outf.split("/")[-1],
        ylabel="TFBS density\nmedian"
    )

    plt.savefig(outf, bbox_inches='tight')


def prep_2x2(tf, arch1, arch2, df):

    # only evaluate enhancers that overlap TF.
    # This excludes enhancers with zero overlaps from background set.
    df = df.loc[df.tfoverlap_bin > 0]

    # split dataframe by two architectures to compare
    dfarch = df.loc[df.arch == arch1]

    if arch2 == "bkgd":
        arch2 = "all_enh_bkgd"
        dfbkgd = df.loc[df.arch != arch1]

    else:
        dfbkgd = df.loc[df.arch == arch2]

    comparison_name = tf + "-" + arch1 + "_v_" + arch2

    # count how many TF overlaps are in each arch.
    TF_in_arch = dfarch.loc[dfarch.tf == tf].shape[0]
    not_TF_in_arch = dfarch.loc[dfarch.tf != tf].shape[0]

    TF_bkgd = dfbkgd.loc[dfbkgd.tf == tf].shape[0]
    not_TF_in_bkgd = dfbkgd.loc[dfbkgd.tf != tf].shape[0]

    a, b, c, d = TF_in_arch, not_TF_in_arch, TF_bkgd, not_TF_in_bkgd

    obs = np.array([[a, b], [c, d]])
    if a+b == 0:
        obs = obs + 1  # add pseudo count

    return obs, comparison_name


def quantify_2x2(obs, comparison_name, min_instances):

    if obs[0][0] > min_instances or obs[1][0] > min_instances:

        OR, P = stats.fisher_exact(obs)
        table = sm.stats.Table2x2(obs)  # get confidence interval
        odds_ci = table.oddsratio_confint()
        newdf = pd.DataFrame({"comparison_name": comparison_name,
                              "a": obs[0][0], "b": obs[0][1],
                              "c": obs[1][0], "d": obs[1][1],
                              "OR": [OR], "P": [P],
                              "ci_lower": [odds_ci[0]],
                              "ci_upper": [odds_ci[1]],
                              })
    else:
        newdf = pd.DataFrame()  # return empty dataframe

    return newdf


def fdr_correction(collection_dict, alpha):

    df = pd.concat(collection_dict.values())

    pvals = df["P"]

    df["reject_null"], df["FDR_P"] = statsmodels.stats.multitest.fdrcorrection(
        pvals, alpha=alpha)

    # other dataframe formatting
    df["arch"] = df["comparison_name"].apply(lambda x: x.split("-")[1])
    df["tf"] = df["comparison_name"].apply(lambda x: x.split("-")[0])
    df["log2"] = np.log2(df["OR"])

    return df


def plot_bar_tf_enrichment(df, cell_line, outf, alpha, taxon2):

    fig, ax = plt.subplots(figsize=(6, 12))
    sns.set("poster")

    x = "tf"
    y = "log2"
    hue = "arch"
    data = df.sort_values(by=y)

    sns.barplot(x=y, y=x, data=data, hue=hue, palette=DERPAL)

    ax.legend(bbox_to_anchor=(1, 1))

    if taxon2 is not None:
        label = cell_line + "_" + taxon2
    else:
        label = cell_line

    ax.set(xlabel=f"OR log2-scale\n FDR<{str(alpha)}", title=label)

    plt.savefig(outf, bbox_inches="tight")


def run_2x2(arch1, arch2, df, min_instances, alpha, taxon2):

    collection_dict = {}

    for tf in df.tf.unique():

        if tf != ".":

            # get 2x2 counts for that TF in arch1, arch2
            obs, comparison_name = prep_2x2(tf, arch1, arch2, df)
            print(obs, comparison_name, tf)

            results = quantify_2x2(obs, comparison_name, min_instances)

            if results.empty is False:
                collection_dict[comparison_name] = results

    # FDR correction
    if len(collection_dict) > 0:  # if there are any results

        results_df = fdr_correction(collection_dict, alpha)

        df = results_df.loc[results_df.reject_null is True]

    else:
        print("\nno results for comparison",
              arch1, "v.", arch2, "in", taxon2)


def run_analysis(cell_line, val, fantombase, encodepath, min_instances, alpha):

    print(cell_line, val)

    fantom, encode, intersection = get_paths(
        cell_line, val, fantombase, encodepath)
    print(fantom, encode, intersection)

    #Bed command
    bed_intersect(fantom, encode, intersection)

    #dataframe
    df = format_df(intersection)

    #get some basic info about Fantom enhancer overlap
    stats_f = os.path.join(RE_DATA, f"stats_{cell_line}.txt") # output file
    archs = ["enh", "syn"]

    for arch in archs:
        count_enhancers(df, arch, stats_f)

        # calculate enhancer TF density
        tf_density, zero = calculate_tf_density(arch, df)

        # non-zero enhancer TF density
        non_zero_tf_density = tf_density.loc[tf_density.tfoverlap_bin > 0]

        # calculate frequency of derived sequences that do not overlap TFBS
        calculate_zero_syn_freq(zero, df, stats_f)

        # plot tf density data
        x, y = "arch", "tf_density"

        # get order of data depending on arch
        order_dict = {"enh":["simple", "complex"],
        "syn":["simple", "complex_core", "complex_derived"]}

        order = order_dict[arch]

        # calculate density with and without zeros
        zeros = {"w_zeros":tf_density, "nonzeros":non_zero_tf_density}

        for kw, data in zeros.items():

            outf = os.path.join(RE, f"{cell_line}_x_encode3_tf_den_{arch}_{kw}.pdf")
            listP, median = mwu(data, arch, kw, stats_f)

            # make a string of the p values for density plotting
            if len(listP) == 2:
                pStr = f"simple_v_core p = {listP[0]},  core_v_der= {listP[1]}"
            else:
                pStr = listP[0]

            plot_bar_tf_density(x, y, data, outf, order, pStr, median)

    ### TF ENRICHMENT IN ARCHITECTURE ###
    # DER V. CORE
    # DER V. BKGD

    # calculate TF enrichment in architecture/syn blocks
    arch1, arch2 = "complex_derived", "complex_core"
    der_v_core = run_2x2(arch1, arch2, df, MIN_INSTANCES, ALPHA, None)

    arch1, arch2 = "complex_derived", "bkgd"
    der_v_bkgd = run_2x2(arch1, arch2, df, MIN_INSTANCES, ALPHA, None)

    arch1, arch2 = "simple", "complex_core"
    simple_v_core = run_2x2(arch1, arch2, df, MIN_INSTANCES, ALPHA, None)

    arch1, arch2 = "simple", "bkgd"
    simple_v_bkgd = run_2x2(arch1, arch2, df, MIN_INSTANCES, ALPHA, None)

    arch1, arch2 = "simple", "complex_derived"
    simple_v_der = run_2x2(arch1, arch2, df, MIN_INSTANCES, ALPHA, None)

    arch1, arch2 = "complex_core", "bkgd"
    core_v_bkgd = run_2x2(arch1, arch2, df, MIN_INSTANCES, ALPHA, None)

    return der_v_core, der_v_bkgd, simple_v_core, simple_v_bkgd, simple_v_der, core_v_bkgd, df


def just_get_df(cell_line, val, fantombase, encodepath,):
    print(cell_line, val)
    fantom, encode, intersection = get_paths(
        cell_line, val, fantombase, encodepath)

    #Bed command
    bed_intersect(fantom, encode, intersection)

    #dataframe
    df = format_df(intersection)

    return df


# %%
sample_dict = get_cell_lines()

# %%
ALPHA = 0.05
MIN_INSTANCES = 10

# %%
cell_line = "HepG2"
val = "ELS_combined_HepG2"

df_file = f"{RE_DATA}{cell_line}_df.tsv"

comp_list = ["der_v_core", "der_v_bkgd",
             "simple_v_core", "simple_v_bkgd", "simple_v_der",
             "core_v_bkgd", "df"]


if os.path.exists(df_file) is False:

    der_v_core, der_v_bkgd,  simple_v_core,\
     simple_v_bkgd, simple_v_der, core_v_bkgd, df = run_analysis(
        cell_line, val, ENHBASE, ENCODEPATH, MIN_INSTANCES, ALPHA)

    data_dict = {  # dictionary of the dataframes and odds ratio comparisons
        "der_v_core": der_v_core,
        "der_v_bkgd": der_v_bkgd,
        "simple_v_core": simple_v_core,
        "simple_v_bkgd": simple_v_bkgd,
        "simple_v_der": simple_v_der,
        "core_v_bkgd": core_v_bkgd,
        "df": df
    }
#%%

#%%
for k, i in data_dict.items():
    print(k, type(i))
    #if bool(i) is True:
    #    print(i.head())
#%%
if bool(data_dict) is True:
    for comp, dataframe in data_dict.items():
        if bool(dataframe) is True:
            outf = f"{RE_DATA}{cell_line}_{comp}.tsv"
            dataframe.to_csv(outf, sep='\t', index=False)

else:
    data_dict = {}
    for comp in comp_list:
        outf = f"{RE_DATA}{cell_line}_{comp}.tsv"
        df = pd.read_csv(outf, sep='\t')
        data_dict[comp] = df

#%%

# %%
df = merge_syn_annot(df)

"""
# IF YOU WANT TO LABEL DERIVED REGIONS
# BY THEIR CORE AGES
# INSTEAD OF THEIR SYNTENIC AGES.
# do the dance - add mrca_2 column, get mrca_2 core age, drop mrca_2 column,
# then add it back,
# but this time to reflect the core_age and core taxon,
# instead of the syntenic age.
"""

REASSIGN_DER_W_CORE_AGE = False


if REASSIGN_DER_W_CORE_AGE is True:
    df = get_core_age(df)
    df = df.drop(["mrca_2"], axis=1)
    df = pd.merge(df, syn[["mrca_2", "taxon2"]], how="left",
                  left_on="core_mrca_2", right_on="mrca_2")
else:
    df = pd.merge(df, syn[["mrca_2", "taxon2"]], how="left")

"""
# IF YOU WANT TO REASSIGN TINY TAXONS
# EUAR -> BORE
# SARG -> VERT
"""

REASSIGN_SMALL_TAXONS = False


if REASSIGN_SMALL_TAXONS is True:
    reassign_small_taxons(df)

# %%
"""
# calculate TF enrichment per age.
"""

mrca_dict = {}

taxon_query = glob.glob(f"{RE_DATA}{cell_line}_*OR_per_MRCA.tsv")

# have you run the taxon enrichments already?
if len(taxon_query) < 6:

    for TAXON2 in df.taxon2.unique():

        print(TAXON2)
        age = df.loc[df.taxon2 == TAXON2].drop_duplicates()

        arch1, arch2 = "complex_derived", "complex_core"
        der_v_core = run_2x2(arch1, arch2, age, MIN_INSTANCES, ALPHA, TAXON2)

        arch1, arch2 = "simple", "complex_core"
        simple_v_core = run_2x2(
            arch1, arch2, age, MIN_INSTANCES, ALPHA, TAXON2)

        arch1, arch2 = "simple", "bkgd"
        simple_v_bkgd = run_2x2(
            arch1, arch2, age, MIN_INSTANCES, ALPHA, TAXON2)

        arch1, arch2 = "complex_derived", "bkgd"
        der_v_bkgd = run_2x2(arch1, arch2, age, MIN_INSTANCES, ALPHA, TAXON2)

        arch1, arch2 = "complex_core", "bkgd"
        core_v_bkgd = run_2x2(arch1, arch2, age, MIN_INSTANCES, ALPHA, TAXON2)

        arch1, arch2 = "complex_derived", "simple"
        der_v_simple = run_2x2(arch1, arch2, age, MIN_INSTANCES, ALPHA, TAXON2)

        results = [der_v_core, simple_v_core, simple_v_bkgd,
                   der_v_bkgd, core_v_bkgd, der_v_simple]

        mrca_dict[TAXON2] = results

# %%


def get_cross_mrca_enrichment(comp, i, mrca_dict):
    comp_dict = {}  # collect TF enrichment across ages
    for key, value in mrca_dict.items():
        compdf = value[i]

        if compdf is None:
            # examples of this is when there are no significant results for an age, or no pairs of complex core and derived (like in vertbrates)
            print("nothing for", key)
        else:
            compdf["taxon2"] = key  # add the taxon annotation to dataframe

            comp_dict[key] = compdf  # core_v_der enrichment df
    # compile the comparison results across ages.
    test = pd.concat(comp_dict.values())  # concat the results
    test = pd.merge(test, syn[["mrca_2", "taxon2"]],
                    how="left")  # merge in MRCA_2 info

    return test


def plot_heatmap(comp, test):
    test = test.drop_duplicates()
    # pivot the results into a table
    table = pd.pivot(test.sort_values(by="mrca_2"),
                     index="tf", columns="mrca_2", values='log2')
    #table = table.dropna(thresh=2)  # drop any Na's
    table = table.replace(-np.Inf, -2)
    table = table.replace(np.Inf, 2)

    if len(table) < 25:
        figsize = (5, 10)
    else:
        figsize = (5, 35)
    # plot
    sns.set("notebook")
    cm = sns.clustermap(table.fillna(0),
                        #mask=(table == 0),
                        cmap="RdBu_r",
                        center=0,
                        robust=True,
                        col_cluster=False,
                        figsize=figsize)

    cm.fig.suptitle(comp)
    outf = f"{RE}{val}_{comp}_clustermap_core_mrca2.pdf"
    plt.savefig(outf, bbox_inches="tight", dpi=300)


#%%
comparison_order = ["der_v_core", "simple_v_core", "simple_v_bkgd",
                    "der_v_bkgd", "core_v_bkgd", "der_v_simple"]

# enumerate and plot results
for i, comp in enumerate(comparison_order):
    print(comp)
    outf = f"{RE_DATA}{cell_line}_{comp}OR_per_MRCA.tsv"
    if os.path.exists(outf) is False:

        test = get_cross_mrca_enrichment(comp, i, mrca_dict)
        test = test.drop_duplicates()
        test["tf"] = test.comparison_name.apply(lambda x: x.split("-")[0])
        test["log2"] = np.log2(test.OR)

        test.to_csv(outf, sep='\t', index=None)
    else:
        test = pd.read_csv(outf, sep='\t')
    plot_heatmap(comp, test)


#%%
synden = data_dict["tf_density_syn"]
synden.head()
test = data_dict["df"]
# add core mrca age
test = pd.merge(test, syn[["mrca", "mrca_2"]])
test = get_core_age(test)
# add age info to tf density dataframe
synden = pd.merge(
    synden, test[["syn_id", "core_mrca_2"]], left_on='id', right_on="syn_id")

#%% plot
x, y = "core_mrca_2", 'tf_density'
hue = "arch"
hue_order = ["simple", "complex_core", "complex_derived"]
xlabs = ["Homo", "Prim", "Euar", "Bore", "Euth",
         "Ther", "Mam", "Amni", "Tetr", "Sarg", "Vert"]
fig, ax = plt.subplots(figsize=(6, 6))
sns.barplot(x, y, data=synden, hue=hue, palette=PAL, hue_order=hue_order)
ax.legend(bbox_to_anchor=(1, 1))
ax.set_xticklabels(xlabs, rotation=90)

outf = f"{RE}{cell_line}_tfbs_density_core_mrca_2.pdf"
outf
plt.savefig(outf, bbox_inches="tight")
#%%

dc = data_dict["der_v_core"]
dc["dummy"] = 0
dc.loc[dc.tf == "FOXA1"]
#%%
dcp = pd.pivot(index='tf', values="log2", columns="dummy", data=dc)
dcp_p = pd.pivot(index='tf', values="reject_null", columns="dummy", data=dc)

#%%
sorted_index = dcp.sort_values(by="dummy", ascending=False).index
dcp_p = dcp_p.reindex(sorted_index)

"""
Reset index to split dataframe
split dataframes into core and derived
reindex core and derived significant arrays by new dataframes
reset TF as the index
"""

d = dcp.reset_index()  # reindex
dp = dcp_p.reset_index()

# split into core and der
core_only, der_only = d.loc[d[0] < 0], d.loc[d[0] >= 0]

# sort values from highest to lowest OR
core_only, der_only = core_only.sort_values(
    by=0), der_only.sort_values(by=0, ascending=False)

# get sorted core index, reindex significance array
coresorted_index = core_only.sort_values(by=0)["tf"].index
dcp_pcore = dp.reindex(coresorted_index)

# get sorted derived index, reindex significance array
dsorted_index = der_only.sort_values(by=0, ascending=False)["tf"].index
dcp_pder = dp.reindex(dsorted_index)

# reset index as TF for all dataframes
core_only = core_only.set_index("tf")
dcp_pcore = dcp_pcore.set_index("tf")

der_only = der_only.set_index("tf")
dcp_pder = dcp_pder.set_index("tf")
#%%

fig, (ax, ax2) = plt.subplots(ncols=2, figsize=(8, 20))
cm = sns.heatmap(core_only,
                 cmap="RdBu_r",
                 center=0,
                 square=True,
                 annot=dcp_pcore.replace({False: "", True: "*"}),
                 fmt="",
                 robust=True,
                 ax=ax
                 )
ax.set(xlabel="", ylabel="", title="core")
cm.set_xlabel(""), cm.set_xticklabels("")

ax2 = sns.heatmap(der_only,
                  cmap="RdBu_r",
                  center=0,
                  square=True,
                  annot=dcp_pder.replace({False: "", True: "*"}),
                  fmt="",
                  robust=True,
                  ax=ax2)
ax2.set(xlabel="", ylabel="", title="derived")
ax2.set_xticklabels("")
outf = os.path.join(RE, "Figure2C_core_derived.pdf")
plt.savefig(outf, bbox_inches="tight")
