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

RE = "/dors/capra_lab/projects/enhancer_ages/landscape/results/cCRE_x_tfbs_encode3/HepG2/pdf/"
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


def just_get_df(cell_line, val, fantombase, encodepath,):
    print(cell_line, val)
    fantom, encode, intersection = get_paths(cell_line, val, fantombase, encodepath)

    #Bed command
    bed_intersect(fantom, encode, intersection)

    #dataframe
    df = format_df(intersection)

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


def add_syn_age_annotation(syn_to_merge):
    SYN_GROUP = "/dors/capra_lab/projects/enhancer_ages/hg38_syn_taxon.bed"
    syn = pd.read_csv(SYN_GROUP, sep = '\t')

    # round all values

    syn[["mrca", "mrca_2"]] = syn[["mrca", "mrca_2"]].round(3)
    syn_to_merge.mrca = syn_to_merge.mrca.round(3)

# do the dance - add mrca_2 column, get mrca_2 core age, drop mrca_2 column, then add it back, but this time to reflect the core_age and core taxon, instead of the syntenic age.
    syn_merged = pd.merge(syn_to_merge, syn[["mrca", "mrca_2"]], how = "left")

    return syn_merged, syn


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


def lump_taxons(df):
    # lump small samples with larger, older ancestors
    df.loc[df.taxon2 == "Sarcopterygian", "taxon2"] = "Vertebrata"
    df.loc[df.taxon2 == "Tetrapoda", "taxon2"] ="Vertebrata"
    df.loc[df.taxon2 == "Euarchontoglires", "taxon2"] = "Boreoeutheria"

    return df

#%% make a dictionary of the cell lines
sample_dict = get_cell_lines()


#%% Load the dataframe
cell_line = "HepG2"
val = "ELS_combined_HepG2"
df = just_get_df(cell_line, val, ENHBASE, ENCODEPATH)
df.head()

#%% how many complex enhancers are there in this dataset?

len(df.loc[(df.core_remodeling ==1) &(df.core ==1), "enh_id"].unique()) #27789

#%%
#get some basic info about enhancer landscapes (simple v. complex), overlap with TFBS

arch = "enh"
totaln, simplen, complexn = count_enhancers(df, arch)

# calculate enhancer TF density
tf_density_enh, zero_enh = calculate_tf_density(arch, df)
zero_enh.id.sum()

# get some basic infor about syntenic landscapes (simple v. core v. derived),
arch = "syn"
totaln, coren, derivedn, simplen = count_enhancers(df, arch)
tf_density_syn, zero_syn = calculate_tf_density(arch, df)
calculate_zero_syn_freq(zero_syn, df, cell_line, RE)


tf_density_syn.head()

#%% evaluate the non-zero syntenic blocks only.

non_zero_syn_tf_density = tf_density_syn.loc[tf_density_syn.tfoverlap_bin>0]
non_zero_syn_tf_density.groupby("arch")["tf_density"].median()

testcore, pcore, test_der, p_der, median = mwu(non_zero_syn_tf_density, arch)


#%% Let's group the dataframe to get syntenic blocks and their TFBS overlap.

syn_ages = df.groupby(["enh_id","syn_id", "syn_len", "arch", "mrca"])["tfoverlap_bin"].sum().reset_index()
syn_ages.head()


#%% #

syn_ages, syn = add_syn_age_annotation(syn_ages) # add MRCA_2 annotations

syn_ages = get_core_age(syn_ages) # annotate core ages for each syntenic block.

syn_ages = syn_ages.rename(columns = {"mrca_2":"syn_mrca_2"})# rename column, preserve the mrca_2 annotation per syntenic block

# add core_mrca_2 taxon annotations
syn_ages = pd.merge(syn_ages, syn[["mrca_2", "taxon2"]], how = "left",
left_on = "core_mrca_2", right_on = "mrca_2")

# include only the values greater than 5
syn_ages = syn_ages.loc[syn_ages.syn_len >5]

# make a boolean for TFBS overlapping syntenic block
syn_ages["tfbs_bool"] = False
syn_ages.loc[syn_ages.tfoverlap_bin >0, "tfbs_bool"] = True
#%% # evaluate zeros as a fraction of the total architecture


zeros_only = syn_ages.loc[syn_ages.tfoverlap_bin ==0]

gz = zeros_only.groupby(["syn_mrca_2", "arch"])["enh_id"].count().reset_index()

gz.columns = ["syn_mrca_2", "arch", 'mrca_zero_counts']

totals = syn_ages.groupby(["arch"])["enh_id"].count().reset_index()
totals.columns = ["arch", "total_arch"]

gz = pd.merge(gz, totals, how = "left")
gz["frac_of_arch"] = gz.mrca_zero_counts.divide(gz.total_arch)
gz.head()


#%% plot
xlabs = ["Homo", "Prim", "Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Sarg", "Vert"]
x, y ="syn_mrca_2", "frac_of_arch"
hue_order = ["simple", "complex_core", "complex_derived"]
fig, ax = plt.subplots()
sns.barplot( data = gz, x= x, y=y, hue = "arch",
hue_order = hue_order, palette = PAL)
ax.set_xticklabels(xlabs, rotation = 90)
ax.set(xlabel = "sequence age",
ylabel = "fraction of total arch")


#%%evaluate zeros as a fraction of the architecture per age

totals_mrcas = syn_ages.groupby(["syn_mrca_2", "arch"])["enh_id"].count().reset_index()
totals_mrcas.columns = ["syn_mrca_2", "arch", "mrca_counts"]

totals_mrcas = pd.merge(totals_mrcas, gz)

totals_mrcas["frac_zero_mrca"] = totals_mrcas.mrca_zero_counts.divide(totals_mrcas.mrca_counts)
totals_mrcas

#%%
xlabs = ["Homo", "Prim", "Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Sarg", "Vert"]

hue_order = ["simple", "complex_core", "complex_derived"]

x, y ="syn_mrca_2", "frac_zero_mrca"
fig, ax = plt.subplots()
sns.barplot( data = totals_mrcas, x= x, y=y, hue = "arch",
hue_order = hue_order, palette = PAL)
ax.set_xticklabels(xlabs, rotation = 90)
ax.set(xlabel = "sequence age",
ylabel = "fraction of arch in mrca")
ax.legend(bbox_to_anchor = (1,1))
outf = f"{RE}zero_frac_per_mrca.pdf"
plt.savefig(outf, bbox_inches = "tight")

#%%
alpha = 0.05

mrca_dict = {}

for arch in syn_ages.arch.unique():
    for mrca_2 in syn_ages.syn_mrca_2.unique():

        in_age = syn_ages.loc[(syn_ages.syn_mrca_2 == mrca_2) & (syn_ages.arch == arch)] # subset to the per age df
        out_age = syn_ages.loc[(syn_ages.syn_mrca_2 != mrca_2) & (syn_ages.arch == arch)] # subset to the per age df
        #print(arch, mrca_2)
        if arch == "complex_derived" and mrca_2 ==0.867:
            continue
        elif arch == "complex_core" and mrca_2 ==0.0:
            continue
        elif arch == "simple" and mrca_2 ==0.0:
            continue
        else:
            ab = in_age.groupby("tfbs_bool")["syn_id"].count().reset_index()
            b,a = ab.iloc[0,1], ab.iloc[1,1]
            cd = out_age.groupby("tfbs_bool")["syn_id"].count().reset_index()
            d,c = cd.iloc[0,1], cd.iloc[1,1]

            obs = [[a,b], [c,d]]

            min_instances = 1
            comparison_name = f"{arch}-{mrca_2}"
            new_df = quantify_2x2(obs, comparison_name, min_instances)
            new_df["mrca_2"] = mrca_2
            new_df["arch"] = arch
            mrca_dict[comparison_name] = new_df


fet_zero_ages = fdr_correction(mrca_dict, alpha) # FDR correction
fet_zero_ages["arch"] = fet_zero_ages.comparison_name.apply(lambda x: x.split("-")[0])
#%%
fet_zero_ages.sort_values(by = "mrca_2")
der_xlabs = xlabs[:-1]
#%%
fet_zero_ages["yerr"]=(fet_zero_ages["ci_upper"] - fet_zero_ages["ci_lower"])
yerrs = fet_zero_ages.pivot(index = "mrca_2", columns = "arch", values = "yerr").fillna(0)

yerrs = yerrs[hue_order]

piv = fet_zero_ages.pivot(index='mrca_2', columns='arch', values='log2')
piv = piv[hue_order].fillna(0)
yerrs
#%%

ax = piv.plot(kind='bar',
yerr=yerrs,
figsize = (6,6), #colormap = PAL,
#ylabel = "Archs that bind TFs enrichment"
)
ax.set_xticklabels(xlabs)
ax.set(
xlabel = "sequence age",
ylabel = "Archs that bind TFs \nOR enrichment")
outf = f"{RE}FigS5_TF_binding_per_arch_age_fet.pdf"
plt.savefig(outf, bbox_inches = 'tight')
#%%
fig, ax = plt.subplots(figsize = (6,6))
x, y = "mrca_2","log2"
data = fet_zero_ages
hue = "arch"
sns.barplot(x = x, y = y, data = data,
hue = hue, palette = PAL, hue_order = hue_order,)
             #yerr= yl)

ax.set(xlabel = "sequence age",
ylabel = "OR architecture binds TF (log2-scaled)",
title = "TFBS enrichment per age, per arch")
ax.set_xticklabels(xlabs)
ax.legend(bbox_to_anchor = (1,1))
outf = f"{RE}FigS5_TF_binding_per_arch_age_fet.pdf"
plt.savefig(outf, bbox_inches = 'tight')
#%%
alpha = 0.05
syn_ages.head()
mrca_dict = {}


    for mrca_2 in syn_ages.syn_mrca_2.unique():

        in_age = syn_ages.loc[(syn_ages.core_mrca_2 == mrca_2) ] # subset to the per age df
        out_age = syn_ages.loc[(syn_ages.core_mrca_2 != mrca_2)] # subset to the per age df
        #print(arch, mrca_2)

        ab = in_age.groupby("tfbs_bool")["syn_id"].count().reset_index()
        b,a = ab.iloc[0,1], ab.iloc[1,1]
        cd = out_age.groupby("tfbs_bool")["syn_id"].count().reset_index()
        d,c = cd.iloc[0,1], cd.iloc[1,1]

        obs = [[a,b], [c,d]]

        min_instances = 1
        comparison_name = f"{arch}-{mrca_2}"
        new_df = quantify_2x2(obs, comparison_name, min_instances)
        new_df["mrca_2"] = mrca_2
        new_df["arch"] = arch
        mrca_dict[comparison_name] = new_df


fet_zero_ages = fdr_correction(mrca_dict, alpha) # FDR correction
fet_zero_ages["arch"] = fet_zero_ages.comparison_name.apply(lambda x: x.split("-")[0])
#%%
fet_zero_ages.sort_values(by = "mrca_2")
der_xlabs = xlabs[:-1]
fig, ax = plt.subplots(figsize = (6,6))
x, y = "mrca_2","log2"
data = fet_zero_ages
#hue = "arch"
sns.barplot(x = x, y = y, data = data, )
#hue = hue, palette = PAL, hue_order = hue_order )
ax.set(xlabel = "sequence age",
ylabel = "OR enhancer binds TF (log2-scaled)",
title = "TFBS enrichment per age")
ax.set_xticklabels(xlabs)
ax.legend(bbox_to_anchor = (1,1))
outf = f"{RE}TF_binding_per_age_fet.pdf"
