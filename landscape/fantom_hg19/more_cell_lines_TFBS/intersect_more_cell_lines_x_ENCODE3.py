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
        ENCODEFILE = cell_line.split("_CL")[0] + ".bed"
    elif cell_line == "all_fantom_enh":
        ENCODEFILE = "trimmed_encRegTfbsClusteredWithCells.liftOver.to.hg19.bed"
    else:
        ENCODEFILE = cell_line + ".bed"

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
        print("previously done enh x encode intersection")


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
        test, p = stats.mannwhitneyu(simple_tfden, complex_tfden)
        print("\n", "simple v. complex enh MWU stat =", round(test,3), "p =", p )

        return test, p, median

    elif arch == "syn":

        simple_tfden = tf_density.loc[tf_density.arch == "simple", "tf_density"]
        core_tfden = tf_density.loc[tf_density.arch == "complex_core", "tf_density"]
        derived_tfden = tf_density.loc[tf_density.arch == "complex_derived", "tf_density"]

        testcore, pcore = stats.mannwhitneyu(simple_tfden, core_tfden)
        print("\n", "simple v. complexcore MWU stat =", round(testcore,3), "p =", pcore)

        test, p = stats.mannwhitneyu(derived_tfden, core_tfden)
        print("\n", "core v. derived MWU stat =", round(test,3), "p =", p)

        return testcore, pcore, test, p, median


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
    zero_overlap = tf_density.loc[tf_density.tfoverlap_bin == 0].groupby("arch")["id"].count()

    return tf_density, zero_overlap


def plot_bar_tf_density(x, y, data, outf, order, p, med):

    # plot enhancer TF density
    fig, ax = plt.subplots(figsize = (6,6))
    sns.set("poster")

    sns.barplot(x, y, data = data, estimator = np.median, order = order)

    ax.set(
    xlabel= "p = %s\n%s" % (p,med),
    title = outf.split("/")[-1],
    ylabel = "TFBS density\nmedian"
    )

    plt.savefig(outf, bbox_inches = 'tight')
    plt.close()


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

        if P<0.05:
            print("\n", comparison_name, obs, round(OR, 2), round(P, 4))

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

    fig, ax = plt.subplots(figsize = (6,9))
    sns.set("poster")

    x = "tf"
    y = "log2"
    hue = "arch"
    data = df

    sns.barplot(x=y, y=x, data=data , hue = hue)

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

    test, p, median = mwu(tf_density_enh, arch)
    plot_bar_tf_density(x, y, data, outf, order, p, median)


    print("\nNon-zero TFBS densities only")
    # plot all enhancer-level data without zeros
    non_zero_tf_density = tf_density_enh.loc[tf_density_enh.tfoverlap_bin>0]

    data = non_zero_tf_density
    outf = make_pdf("%s_enh_x_encode3_tf_density_%s_non_zero_tf_density" % (cell_line, arch), RE)
    test, p, median = mwu(non_zero_tf_density, arch)
    plot_bar_tf_density(x, y, data, outf, order, p, median)



    # calculate syn-level TF density
    print("\nSyn TFBS densities")
    arch = "syn"
    totaln, coren, derivedn, simplen = count_enhancers(df, arch)
    tf_density_syn, zero_syn = calculate_tf_density(arch, df)


    # plot syn block TF density
    order = ["simple", "complex_core", "complex_derived"]
    data = tf_density_syn
    outf = make_pdf("%s_enh_x_encode3_tf_density_%s" % (cell_line, arch), RE)

    testcore, pcore, test, p, median = mwu(tf_density_syn, arch)
    new_p = "core_v_der p = %s, simple_v_core = %s" %(pcore, p)
    plot_bar_tf_density(x, y, data, outf, order, new_p, median)


    print("\nNon-zero syn TFBS densities only")
    non_zero_syn_tf_density = tf_density_syn.loc[tf_density_syn.tfoverlap_bin>0]

    # plot syn block TF density
    data = non_zero_syn_tf_density
    outf = make_pdf("%s_syn_x_encode3_tf_density_%s_non_zero_tf_density" % (cell_line, arch), RE)

    testcore, pcore, test, p, median = mwu(non_zero_syn_tf_density, arch)
    new_p = "core_v_der p = %s, simple_v_core = %s" %(pcore, p)
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

#%%

cell_line = "all_fantom_enh"
val = sample_dict[cell_line]


der_v_core, der_v_bkgd, tf_density_enh, tf_density_syn, df = run_analysis(cell_line, val, FANTOMBASE, ENCODEPATH, MIN_INSTANCES, ALPHA)

results_dict[cell_line] = df
tf_den_enh[cell_line] = tf_density_enh
tf_den_syn[cell_line] = tf_density_syn
der_v_bkgd_dict[cell_line] = der_v_bkgd
der_v_core_dict[cell_line] = der_v_core


#%%
ALPHA = 0.1
MIN_INSTANCES = 20
#%%

cell_line = "HepG2"
val = sample_dict[cell_line]


der_v_core, der_v_bkgd, tf_density_enh, tf_density_syn, df = run_analysis(cell_line, val, FANTOMBASE, ENCODEPATH, MIN_INSTANCES, ALPHA)

results_dict[cell_line] = df
tf_den_enh[cell_line] = tf_density_enh
tf_den_syn[cell_line] = tf_density_syn
der_v_bkgd_dict[cell_line] = der_v_bkgd
der_v_core_dict[cell_line] = der_v_core


der_v_core.head()


#%% run the CL file
cell_line = "HepG2_CL"
val = sample_dict[cell_line]


der_v_core, der_v_bkgd, tf_density_enh, tf_density_syn, df = run_analysis(cell_line, val, FANTOMBASE, ENCODEPATH, MIN_INSTANCES, ALPHA)

results_dict[cell_line] = df
tf_den_enh[cell_line] = tf_density_enh
tf_den_syn[cell_line] = tf_density_syn
der_v_bkgd_dict[cell_line] = der_v_bkgd
der_v_core_dict[cell_line] = der_v_core


#%%

cell_line = "K562"
val = sample_dict[cell_line]

der_v_core, der_v_bkgd, tf_density_enh, tf_density_syn, df = run_analysis(cell_line, val, FANTOMBASE, ENCODEPATH, MIN_INSTANCES, ALPHA)

results_dict[cell_line] = df
tf_den_enh[cell_line] = tf_density_enh
tf_den_syn[cell_line] = tf_density_syn
der_v_bkgd_dict[cell_line] = der_v_bkgd
der_v_core_dict[cell_line] = der_v_core

#%%

cell_line = "GM12878"
val = sample_dict[cell_line]

der_v_core, der_v_bkgd, tf_density_enh, tf_density_syn, df = run_analysis(cell_line, val, FANTOMBASE, ENCODEPATH, MIN_INSTANCES, ALPHA)

results_dict[cell_line] = df
tf_den_enh[cell_line] = tf_density_enh
tf_den_syn[cell_line] = tf_density_syn
der_v_bkgd_dict[cell_line] = der_v_bkgd
der_v_core_dict[cell_line] = der_v_core
#%%
cell_line = "GM12878_CL"
val = sample_dict[cell_line]

der_v_core, der_v_bkgd, tf_density_enh, tf_density_syn, df = run_analysis(cell_line, val, FANTOMBASE, ENCODEPATH, MIN_INSTANCES, ALPHA)

results_dict[cell_line] = df
tf_den_enh[cell_line] = tf_density_enh
tf_den_syn[cell_line] = tf_density_syn
der_v_bkgd_dict[cell_line] = der_v_bkgd
der_v_core_dict[cell_line] = der_v_core


#%%

cell_line = "A549"
val = sample_dict[cell_line]

der_v_core, der_v_bkgd, tf_density_enh, tf_density_syn, df = run_analysis(cell_line, val, FANTOMBASE, ENCODEPATH, MIN_INSTANCES, ALPHA)

results_dict[cell_line] = df
tf_den_enh[cell_line] = tf_density_enh
tf_den_syn[cell_line] = tf_density_syn
der_v_bkgd_dict[cell_line] = der_v_bkgd
der_v_core_dict[cell_line] = der_v_core


#%%

cell_line = "liver"
val = sample_dict[cell_line]

der_v_core, der_v_bkgd, tf_density_enh, tf_density_syn, df = run_analysis(cell_line, val, FANTOMBASE, ENCODEPATH, MIN_INSTANCES, ALPHA)

results_dict[cell_line] = df
tf_den_enh[cell_line] = tf_density_enh
tf_den_syn[cell_line] = tf_density_syn
der_v_bkgd_dict[cell_line] = der_v_bkgd
der_v_core_dict[cell_line] = der_v_core
