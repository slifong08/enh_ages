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

ENCODEPATH = "/dors/capra_lab/data/encode/encode3_hg38/TF/"

RE = "/dors/capra_lab/projects/enhancer_ages/landscape/results/cCRE_x_tfbs_encode3/HepG2/"
RE_DATA = RE + "data/"

if os.path.exists(RE_DATA) == False:
    os.mkdir(RE_DATA)

colors = [ "amber", "dusty purple", "windows blue"]
PAL = sns.xkcd_palette(colors)
sns.palplot(PAL)

colors = [ "windows blue"]
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
    FANTOMFILE = "syn_breaks_%s_ages.bed" % file_tag
    FANTOM = os.path.join(FANTOMPATH, FANTOMFILE)

    if "CL" in cell_line:
        ENCODEFILE = "cells/" + cell_line.split("_CL")[0] + ".bed.gz"
    elif cell_line == "all_fantom_enh":
        ENCODEFILE = "trimmed_encRegTfbsClusteredWithCells.bed"
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
        print("previously done enh x encode intersection")


def get_core_age(df):
    core = df.groupby("enh_id")["mrca_2"].max().reset_index()
    core.columns = ["enh_id", 'core_mrca_2']
    df = pd.merge(df, core, how = "left")

    return df


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


    checkrowone = a + b

    if checkrowone > 0:
        obs = [[a,b], [c,d]]

        return obs, comparison_name
    else:
        print("no obs for", tf)

        obs = [[0,0], [0,0]]

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
    else:
        newdf = pd.DataFrame() # return empty dataframe


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


def plot_bar_tf_enrichment(df, cell_line, outf, alpha, taxon2):

    fig, ax = plt.subplots(figsize = (6,12))
    sns.set("poster")

    x = "tf"
    y = "log2"
    hue = "arch"
    data = df.sort_values(by = y)

    sns.barplot(x=y, y=x, data=data , hue = hue, palette = DERPAL)

    ax.legend(bbox_to_anchor = (1,1))

    if taxon2 !=None:
        label = cell_line + "_" + taxon2
    else:
        label = cell_line

    ax.set(xlabel = "OR log2-scale\n FDR<%s" % str(alpha), title = label)

    plt.savefig(outf, bbox_inches = "tight")


def run_2x2(arch1, arch2, df, min_instances, alpha, taxon2):

    collection_dict = {}

    for tf in df.tf.unique():

        if tf != ".":

            obs, comparison_name = prep_2x2(tf, arch1, arch2, df)

            results = quantify_2x2(obs, comparison_name, min_instances)

            if results.empty ==False:
                collection_dict[comparison_name] = results

    #FDR correction
    if len(collection_dict) > 0: # if there are any results

        results_df = fdr_correction(collection_dict, alpha)

        df = results_df.loc[results_df.reject_null == True]


        if len(df)>0: # if there are any significant results, plot them!

            if taxon2 != None:
                outf = make_pdf("%s_enh_x_encode3_sig_tf_arch_enrichment_%s_v_%s_FDR_%s_%s" % (cell_line, arch1, arch2, alpha, taxon2), RE)
                #plot_bar_tf_enrichment(df, cell_line, outf, alpha, taxon2)

            else:
                outf = make_pdf("%s_enh_x_encode3_sig_tf_arch_enrichment_%s_v_%s_FDR_%s" % (cell_line, arch1, arch2, alpha), RE)
                #plot_bar_tf_enrichment(df, cell_line, outf, alpha, taxon2)

            return results_df

        else:
                print("\nno sig results for comparison", arch1, "v.", arch2, "in", taxon2)

    else:
        print("\nnot any results for comparison", arch1, "v.", arch2, "in", taxon2)


def run_analysis(cell_line, val, fantombase, encodepath, min_instances, alpha):

    print(cell_line, val)
    fantom, encode, intersection = get_paths(cell_line, val, fantombase, encodepath)
    print(fantom, encode, intersection)
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

    return der_v_core, der_v_bkgd, tf_density_enh, tf_density_syn, simple_v_core, simple_v_bkgd, simple_v_der, core_v_bkgd, df


def just_get_df(cell_line, val, fantombase, encodepath,):
    print(cell_line, val)
    fantom, encode, intersection = get_paths(cell_line, val, fantombase, encodepath)

    #Bed command
    bed_intersect(fantom, encode, intersection)

    #dataframe
    df = format_df(intersection)

    return df

#%%
sample_dict = get_cell_lines()

#%%
ALPHA = 0.1
MIN_INSTANCES =10

#%%
cell_line = "HepG2"
val = "ELS_combined_HepG2"
#val = sample_dict[cell_line]
df_file = f"{RE_DATA}{cell_line}_df.tsv"

comp_list = ["der_v_core", "der_v_bkgd", "tf_density_enh", "tf_density_syn",
"simple_v_core",
"simple_v_bkgd", "simple_v_der", "core_v_bkgd", "df"]


if os.path.exists(df_file) == False:

    der_v_core, der_v_bkgd, tf_density_enh, tf_density_syn, simple_v_core, simple_v_bkgd, simple_v_der, core_v_bkgd, df = run_analysis(cell_line, val, ENHBASE, ENCODEPATH, MIN_INSTANCES, ALPHA)

    data_dict = {  # dictionary of the dataframes and odds ratio comparisons
    "der_v_core":der_v_core,
    "der_v_bkgd":der_v_bkgd,
    "tf_density_enh":tf_density_enh,
    "tf_density_syn":tf_density_syn,
    "simple_v_core":simple_v_core,
    "simple_v_bkgd":simple_v_bkgd,
    "simple_v_der": simple_v_der,
    "core_v_bkgd":core_v_bkgd,
    "df":df
    }

    for comp, dataframe in data_dict.items():
        outf = f"{RE_DATA}{cell_line}_{comp}.tsv"
        dataframe.to_csv(outf, sep = '\t', index = False)
else:
    data_dict = {}
    for comp in comp_list:
        outf = f"{RE_DATA}{cell_line}_{comp}.tsv"
        df = pd.read_csv(outf, sep = '\t')
        data_dict[comp] = df


#%%

SYN_GROUP = "/dors/capra_lab/projects/enhancer_ages/hg38_syn_taxon.bed"
syn = pd.read_csv(SYN_GROUP, sep = '\t')

# round all values


syn[["mrca", "mrca_2"]] = syn[["mrca", "mrca_2"]].round(3)
df.mrca = df.mrca.round(3)

# do the dance - add mrca_2 column, get mrca_2 core age, drop mrca_2 column, then add it back, but this time to reflect the core_age and core taxon, instead of the syntenic age.
df = pd.merge(df, syn[["mrca", "mrca_2"]], how = "left")
df = get_core_age(df)
df = df.drop(["mrca_2"], axis = 1)
df = pd.merge(df, syn[["mrca_2", "taxon2"]], how = "left", left_on = "core_mrca_2", right_on = "mrca_2")
df.head(10)



# lump small samples with larger, older ancestors
df.loc[df.taxon2 == "Sarcopterygian", "taxon2"] = "Vertebrata"
df.loc[df.taxon2 == "Tetrapoda", "taxon2"] ="Vertebrata"
df.loc[df.taxon2 == "Euarchontoglires", "taxon2"] = "Boreoeutheria"

#%%

mrca_dict = {}
# calculate TF enrichment in architecture/syn blocks

for TAXON2 in df.taxon2.unique():

    print(TAXON2)
    age = df.loc[df.taxon2 == TAXON2]

    arch1, arch2 = "complex_derived", "complex_core"
    der_v_core = run_2x2(arch1, arch2, age, MIN_INSTANCES, ALPHA, TAXON2)


    arch1, arch2 = "simple", "complex_core"
    simple_v_core = run_2x2(arch1, arch2, age, MIN_INSTANCES, ALPHA, TAXON2)

    arch1, arch2 = "simple", "bkgd"
    simple_v_bkgd = run_2x2(arch1, arch2, age, MIN_INSTANCES, ALPHA, TAXON2)

    arch1, arch2 = "complex_derived", "bkgd"
    der_v_bkgd = run_2x2(arch1, arch2, age, MIN_INSTANCES, ALPHA, TAXON2)

    arch1, arch2 = "complex_core", "bkgd"
    core_v_bkgd = run_2x2(arch1, arch2, age, MIN_INSTANCES, ALPHA, TAXON2)

    arch1, arch2 = "complex_derived", "simple"
    der_v_simple = run_2x2(arch1, arch2, age, MIN_INSTANCES, ALPHA, TAXON2)

    results = [der_v_core, simple_v_core, simple_v_bkgd, der_v_bkgd, core_v_bkgd, der_v_simple]

    mrca_dict[TAXON2] = results

#%%

def get_cross_mrca_enrichment(comp, i, mrca_dict):
    comp_dict = {} # collect TF enrichment across ages
    for key, value in mrca_dict.items():
        compdf = value[i]

        if compdf is None:
            print("nothing for", key) # examples of this is when there are no significant results for an age, or no pairs of complex core and derived (like in vertbrates)
        else:
            compdf["taxon2"] = key # add the taxon annotation to dataframe

            comp_dict[key] = compdf # core_v_der enrichment df
    # compile the comparison results across ages.
    test = pd.concat(comp_dict.values()) # concat the results
    test = pd.merge(test, syn[["mrca_2", "taxon2"]], how = "left") # merge in MRCA_2 info

    return test

def plot_heatmap(comp, test):
    test = test.drop_duplicates()
    # pivot the results into a table
    table = pd.pivot(test.loc[test.reject_null==True].sort_values(by = "mrca_2"),
    index = "tf", columns = [ "mrca_2", "taxon2"], values = 'log2') # pivot only the significant results
    table = table.dropna(thresh = 0) # drop any Na's
    if len(table)<25:
        figsize = (5,10)
    else:
        figsize= (5,30)
    # plot
    sns.set("notebook")
    cm = sns.clustermap(table.fillna(0), mask = (table==0),
     cmap = "RdBu_r", center = 0, col_cluster = False,
    figsize = figsize)

    cm.fig.suptitle(comp)
    outf = f"{RE}{val}_{comp}_clustermap_core_mrca2.pdf"
    plt.savefig(outf, bbox_inches = "tight", dpi = 300)

#%%
comparison_order=["der_v_core" ,"simple_v_core","simple_v_bkgd",
"der_v_bkgd", "core_v_bkgd", "der_v_simple"]

# enumerate and plot results
for i, comp in enumerate(comparison_order):
    print(comp)
    test = get_cross_mrca_enrichment(comp, i, mrca_dict)
    outf = f"{RE_DATA}{cell_line}_{comp}OR_per_MRCA.tsv"
    test.to_csv(outf, sep = '\t', index = None)
    plot_heatmap(comp, test)

#%%
synden = data_dict["tf_density_syn"]
synden.head()
test = data_dict["df"]
# add core mrca age
test = pd.merge(test, syn[["mrca", "mrca_2"]])
test = get_core_age(test)
# add age info to tf density dataframe
synden = pd.merge(synden, test[["syn_id", "core_mrca_2"]], left_on = 'id', right_on = "syn_id")

#%% plot
x ,y = "core_mrca_2", 'tf_density'
hue = "arch"
hue_order = ["simple", "complex_core", "complex_derived"]
xlabs = ["Homo", "Prim", "Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Sarg", "Vert"]
fig, ax = plt.subplots(figsize = (6,6))
sns.barplot(x, y, data = synden, hue = hue, palette = PAL, hue_order = hue_order)
ax.legend(bbox_to_anchor = (1,1))
ax.set_xticklabels(xlabs, rotation = 90)

outf = f"{RE}{cell_line}_tfbs_density_core_mrca_2.pdf"
outf
plt.savefig(outf, bbox_inches = "tight")
