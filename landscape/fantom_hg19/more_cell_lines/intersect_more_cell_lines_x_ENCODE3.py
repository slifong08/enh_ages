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

RE = "/dors/capra_lab/projects/enhancer_ages/fantom/results/encode3/"



#%% Functions


def get_cell_lines():
    sample_dict = {
    "GM12878": "GM12878_FANTOM5_tpm_hg19",
    "HepG2": "HEPG2_FANTOM5_tpm_hg19",
    "A549":"A549_FANTOM5_tpm_hg19",
    #"H1-hESC":"H1-esc_FANTOM5_tpm_hg19",
    #"PC-3": "PC-3_FANTOM5_tpm_hg19",
    "K562":"CL_0000094_granulocyte_expressed_enhancers",
    #"MCF7":"CL_0002327_mammary_epithelial_cell_expressed_enhancers",
    "liver":"UBERON_0002107_liver_expressed_enhancers",
    }

    return sample_dict


def get_paths(cell_line, file_tag, fantombase, encodepath):

    FANTOMPATH = os.path.join(fantombase, file_tag, "ages")
    FANTOMFILE = "syn_breaks_%s_ages.bed" % file_tag
    FANTOM = os.path.join(FANTOMPATH, FANTOMFILE)

    ENCODEFILE = cell_line + ".bed"
    ENCODE = os.path.join(encodepath, ENCODEFILE)


    INTERSECTIONPATH = os.path.join(fantombase, file_tag)
    INTERSECTIONFILE = "%s_x_ENCODE.bed" % file_tag
    INTERSECTION = os.path.join(INTERSECTIONPATH, INTERSECTIONFILE)

    return FANTOM, ENCODE, INTERSECTION


def bed_intersect(fantom, encode, intersection):

    cmd = "bedtools intersect -a %s -b %s -wao > %s" % (fantom, encode, intersection)

    subprocess.call(cmd, shell = True)

    print(cmd)


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
    header = None)

    df.columns = cols # add column names

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
    mean = (tf_density.groupby("arch")["tf_density"].mean())
    print("means", mean)

    if arch == "enh":

        #stratify dataframe by simple and complex arch
        simple_tfden = tf_density.loc[tf_density.arch == "simple", "tf_density"]
        complex_tfden = tf_density.loc[tf_density.arch == "complex", "tf_density"]

        # calculate MWU
        test, p = stats.mannwhitneyu(simple_tfden, complex_tfden)
        print("simple v. complex enh MWU stat =", round(test,3), "p =", p )

        return test, p, mean

    elif arch == "syn":

        simple_tfden = tf_density.loc[tf_density.arch == "simple", "tf_density"]
        core_tfden = tf_density.loc[tf_density.arch == "complex_core", "tf_density"]
        derived_tfden = tf_density.loc[tf_density.arch == "complex_derived", "tf_density"]

        testcore, pcore = stats.mannwhitneyu(simple_tfden, core_tfden)
        print("simple v. complexcore MWU stat =", round(testcore,3), "p =", pcore)

        test, p = stats.mannwhitneyu(derived_tfden, core_tfden)
        print("core v. derived MWU stat =", round(test,3), "p =", p)

        return testcore, pcore, test, p, mean


def calculate_tf_density(arch):

    density_cols = ["id", "len", "arch", "tfoverlap_bin", "tf_density"]

    if arch == "enh":

        tf_density = df.groupby(["enh_id", "enh_len", "overallarch"])["tfoverlap_bin"].sum().reset_index().drop_duplicates()
        tf_density["tf_density"] = tf_density["tfoverlap_bin"].divide(tf_density.enh_len)


    elif arch == "syn":

        tf_density = df.groupby(["syn_id", "syn_len", "arch"])["tfoverlap_bin"].sum().reset_index()
        tf_density["tf_density"] = tf_density["tfoverlap_bin"].divide(tf_density.syn_len)

    # rename columns
    tf_density.columns = density_cols

    # how many enhancers do not overlap  TFs?
    zero_overlap = tf_density.loc[tf_density.tfoverlap_bin == 0].groupby("arch")["id"].count()

    return tf_density, zero_overlap


def plot_bar_tf_density(x, y, data, outf, order):

    # plot enhancer TF density (including zeros)

    fig, ax = plt.subplots(figsize = (6,6))
    sns.set("poster")

    sns.barplot(x, y, data = data, estimator = np.mean, order = order)
    ax.set_xlabel(outf.split("/")[-1])
    plt.savefig(outf, bbox_inches = 'tight')


def make_pdf(file_name, RE):

    OUTFILE = file_name + ".pdf"
    OUTF = os.path.join(RE, OUTFILE)

    return OUTF


def prep_2x2(tf, arch1, arch2, df):

    # only evaluate enhancers that overlap TF
    df = df.loc[df.tf != "."]

    # split dataframe by architectures to compare
    dfarch = df.loc[df.arch == arch1]
    dfbkgd = df.loc[df.arch == arch2]

    # count how many TF overlaps are in each arch.
    TF_in_arch = len(dfarch.loc[dfarch.tf == tf])
    TF_bkgd = len(dfbkgd.loc[dfbkgd.tf == tf])
    not_TF_in_arch = len(dfarch.loc[dfarch.tf != tf])
    not_TF_in_bkgd = len(dfbkgd.loc[dfbkgd.tf != tf])

    a, b, c, d = TF_in_arch, not_TF_in_arch, TF_bkgd, not_TF_in_bkgd

    obs = [[a,b], [c,d]]

    return obs


def quantify_2x2(obs, comparison_name):

    if obs[0][0] > 100 or obs[1][0]>100:

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

    df["reject_null"], df["FDR_P"] = statsmodels.stats.multitest.fdrcorrection(pvals, alpha=0.1)

    # other dataframe formatting
    df["arch"] = df["comparison_name"].apply(lambda x: x.split("-")[1])
    df["tf"] = df["comparison_name"].apply(lambda x: x.split("-")[0])
    df["log2"]= np.log2(df["OR"])

    return df


def plot_bar_tf_enrichment(sigresults, cell_line, outf):

    fig, ax = plt.subplots(figsize = (6,9))
    sns.set("poster")

    x = "tf"
    y = "log2"
    hue = "arch"
    data = sigresults

    sns.barplot(x=y, y=x, data=data , hue = hue, )

    ax.legend(bbox_to_anchor = (1,1))
    ax.set(xlabel = "OR log2-scale\n FDR<10%", title = cell_line)

    plt.savefig(outf, bbox_inches = "tight")


def run_analysis(cell_line, val, FANTOMBASE, ENCODEPATH):


    print(cell_line, val)
    FANTOM, ENCODE, INTERSECTION = get_paths(cell_line, val, FANTOMBASE, ENCODEPATH)

    #Bed command
    bed_intersect(FANTOM, ENCODE, INTERSECTION)

    #dataframe
    df = format_df(INTERSECTION)

    #get some basic info about Fantom enhancer overlap
    arch = "enh"
    totaln, simplen, complexn = count_enhancers(df, arch)

    # calculate enhancer TF density
    tf_density_enh, zero_enh = calculate_tf_density(arch)

    # plot all enhancer-level data
    x = "arch"
    y = "tf_density"
    order = ["simple", "complex"]
    data = tf_density_enh
    outf = make_pdf("%s_enh_x_encode3_tf_density_%s"  % (cell_line, arch), RE)

    plot_bar_tf_density(x, y, data, outf, order)
    mwu(tf_density_enh, arch)

    # plot all enhancer-level data without zeros
    plot_data = tf_density_enh.loc[tf_density_enh.tfoverlap_bin>0]

    data = plot_data
    outf = make_pdf("%s_enh_x_encode3_tf_density_%s_no_zeros" % (cell_line, arch), RE)
    plot_bar_tf_density(x, y, data, outf, order)
    mwu(plot_data, arch)


    # calculate syn-level TF density
    arch = "syn"
    totaln, coren, derivedn, simplen = count_enhancers(df, arch)
    tf_density_syn, zero_syn = calculate_tf_density(arch)


    # plot syn block TF density
    x = "arch"
    y = "tf_density"
    order = ["simple", "complex_core", "complex_derived"]
    data = tf_density_syn
    outf = make_pdf("%s_enh_x_encode3_tf_density_%s" % (cell_line, arch), RE)

    plot_bar_tf_density(x, y, data, outf, order)
    mwu(tf_density_syn, arch)

    # calculate TF enrichment in architecture/syn blocks
    collection_dict = {}

    for tf in df.tf.unique():
        arch1 = "derived"
        arch2 = "complex_core"
        print(arch1, arch2)
        if tf != ".":
            comparison_name = tf + "-" + arch

            obs = prep_2x2(tf, arch1, arch2, df)

            results = quantify_2x2(obs, comparison_name)

            collection_dict[comparison_name] = results
    results_df = fdr_correction(collection_dict)
    sigresults = results_df.loc[results_df.reject_null == True]
    outf = make_pdf("%s_enh_x_encode3_sig_tf_arch_enrichment_%s" % (cell_line, arch), RE)

    return results_df, tf_density_enh, tf_density_syn, df
    #plot_bar_tf_enrichment(sigresults, cell_line, outf)




#%%
sample_dict = get_cell_lines()


tf_arch_dict = {} # collect all the dataframes for tf enrichment
tf_den_enh = {}
tf_den_syn = {}

#%%


cell_line = "HepG2"
val = sample_dict[cell_line]

sigresults, tf_density_enh, tf_density_syn, df = run_analysis(cell_line, val, FANTOMBASE, ENCODEPATH)

results_dict[cell_line] = sigresults
tf_den_enh[cell_line] = tf_density_enh
tf_den_syn[cell_line] = tf_density_syn


#%%
sigresults.head()
#%%

cell_line = "K562"
val = sample_dict[cell_line]

sigresults, tf_density_enh, tf_density_syn, df = run_analysis(cell_line, val, FANTOMBASE, ENCODEPATH)

results_dict[cell_line] = sigresults
tf_den_enh[cell_line] = tf_density_enh
tf_den_syn[cell_line] = tf_density_syn

#%%

cell_line = "GM12878"
val = sample_dict[cell_line]

sigresults, tf_density_enh, tf_density_syn, df = run_analysis(cell_line, val, FANTOMBASE, ENCODEPATH)

results_dict[cell_line] = sigresults
tf_den_enh[cell_line] = tf_density_enh
tf_den_syn[cell_line] = tf_density_syn
#%%

cell_line = "A549"
val = sample_dict[cell_line]

sigresults, tf_density_enh, tf_density_syn, df = run_analysis(cell_line, val, FANTOMBASE, ENCODEPATH)

results_dict[cell_line] = sigresults
tf_den_enh[cell_line] = tf_density_enh
tf_den_syn[cell_line] = tf_density_syn
#%%

cell_line = "liver"
val = sample_dict[cell_line]

sigresults, tf_density_enh, tf_density_syn, df = run_analysis(cell_line, val, FANTOMBASE, ENCODEPATH)

results_dict[cell_line] = sigresults
tf_den_enh[cell_line] = tf_density_enh
tf_den_syn[cell_line] = tf_density_syn
