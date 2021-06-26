import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy import stats
import seaborn as sns
import statsmodels
import statsmodels.api as sm
import subprocess
#%%

FANTOMPATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/non-genic/ages/"
FANTOMFILE = "syn_breaks_no-exon_all_fantom_enh_ages.bed"
FANTOM = os.path.join(FANTOMPATH, FANTOMFILE)

#FANTOM_TFBS_ONLY = f"{FANTOMPATH}enh_tfbs_only.txt"

SHUFPATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/shuffle/first_round_breaks"
SHUFFILE = "noexon.bed"
SHUFF = os.path.join(SHUFPATH, SHUFFILE)


RE = "/dors/capra_lab/projects/enhancer_ages/landscape/results/fantom/arch_features/all_arch/"


#%%


def add_arch_labels(df):

    df["arch"] = "complex_core"

    df.loc[df.core_remodeling ==0, "arch"] = "simple"
    df.loc[df.core ==0, "arch"] = "complex_derived"

    return df


def format_syndf(enh_age_file):

    syn_cols = ["chr_syn", "start_syn", "end_syn",
    "enh_id",
    "chr", "start", "end",
    "seg_index", "core_remodeling", "core",
    "mrca",]

    syn = pd.read_csv(enh_age_file, sep ='\t', header = None, names = syn_cols)

    syn["syn_id"] = syn.chr_syn + ":" + syn.start_syn.map(str) + "-" + syn.end_syn.map(str)

    syn["syn_len"] = syn.end_syn - syn.start_syn
    syn["enh_len"] = syn.end - syn.start


    # age and taxon file
    syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
    syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
    syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

    syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2"]] # whittle down the df
    syn["mrca"] = syn["mrca"].round(3) # round the ages

    syn = pd.merge(syn, syn_gen_bkgd, how = "left", on = "mrca")

    labeled_syn = add_arch_labels(syn) # add architecture labels

    return labeled_syn


def get_core_der_frac(enh):
    cor_dict = {}
    mrca_lst =  enh.mrca_2.unique()
    for mrca_2 in mrca_lst:
        if mrca_2 != 0.00:
            # get the set of core enhancer ages
            core = set(enh.loc[(enh["core"] ==1)
            & (enh["core_remodeling"] == 1)
            & (enh["mrca_2"] == mrca_2), "enh_id"])

            collection_results = {} # collect all the n
            # per derived age
            for mrca_2_der in mrca_lst:

                # get the enhancer ids with the core age and with derived age
                der = set(enh.loc[(enh["enh_id"].isin(core)) & (enh['mrca_2']== mrca_2_der), "enh_id"])

                # quantify % of cores that have der region of mrca_2_der age
                percent_overlap = len(der)/len(core)

                results = pd.DataFrame({
                "core_mrca":[mrca_2],
                "der_mrca":[mrca_2_der],
                "core_enh_count":[len(core)],
                "der_enh_count":[len(der)],
                "core_der_overlap":[percent_overlap]})

                collection_results[mrca_2_der] = results
            all_results = pd.concat(collection_results.values())

            cor_dict[mrca_2] = all_results


    corr = pd.concat(cor_dict.values())
    return corr


def calc_OR(corr, shuf_corr):

    collection_dict = {}

    # iterate through core and derived ages
    for core_mrca_2 in corr.core_mrca.unique():

        for der_mrca_2 in corr.der_mrca.unique():
            # only test core /derived pairs according to definition that
            # derived sequences are younger than core sequences
            print(core_mrca_2, "-", der_mrca_2)
            if der_mrca_2 < core_mrca_2 and core_mrca_2 != 0.0:

                comparison_name = f"core-{core_mrca_2}_der-{der_mrca_2}"
                total_core = corr.loc[corr["core_mrca"] == core_mrca_2, "core_enh_count"].iloc[0]
                total_core_shuf = shuf_corr.loc[shuf_corr["core_mrca"] == core_mrca_2, "core_enh_count"].iloc[0]

                a = corr.loc[(corr["core_mrca"] == core_mrca_2)&(corr["der_mrca"] == der_mrca_2), "der_enh_count"].iloc[0] # the cores that overlap derived regions of an age
                b = total_core - a # the cores that do not overlap derived regions of an age
                c = shuf_corr.loc[(shuf_corr["core_mrca"] == core_mrca_2)&(shuf_corr["der_mrca"] == der_mrca_2), "der_enh_count"].iloc[0] # the shuffled cores that overlap derived regions of an age
                d = total_core_shuf - c # the shuffled cores that do not overlap derived regions of an age

                obs = [[a,b], [c,d]]


                OR, P = stats.fisher_exact(obs)
                table = sm.stats.Table2x2(obs) # get confidence interval
                odds_ci = table.oddsratio_confint()
                newdf = pd.DataFrame({"comparison_name":comparison_name,
                                      "a":obs[0][0], "b":obs[0][1],
                                      "c":obs[1][0], "d":obs[1][1],
                                      "OR":[OR], "P":[P],
                                      "ci_lower" :[odds_ci[0]],
                                      "ci_upper" :[odds_ci[1]],
                                      "core_mrca": [core_mrca_2],
                                      "der_mrca": [der_mrca_2]
                                    })
                collection_dict[comparison_name] = newdf
                if P<0.05:
                    print(comparison_name, obs, round(OR, 2), round(P, 4))

    # do an FDR correction
    fdr_corrected = fdr_correction(collection_dict)
    # calculate log10p for enrichment/depletion
    fdr_corrected["-log10p"] = np.log10(fdr_corrected["FDR_P"]) *(-1)

    # calculate log2 OR
    fdr_corrected["log2"] = np.log2(fdr_corrected["OR"])


    return fdr_corrected


def fdr_correction(collection_dict):

    df = pd.concat(collection_dict.values())

    pvals = df["P"]

    df["reject_null"], df["FDR_P"] = statsmodels.stats.multitest.fdrcorrection(pvals, alpha=0.05)

    return df


def add_zeros_to_table(table):


    table[0] = 0
    table = table.sort_index()
    if 0.957 not in list(table.index):
        table.loc[0.957] = 0

    table = table.replace(-np.Inf, np.nan) # clean up and fill negative infinitis
    table = table.fillna(0)
    table = table.round(2)
    table = table[[ 0.0, 0.126, 0.131, 0.152, 0.175, 0.308, 0.38, 0.49, 0.656, 0.957,]]

    return table


def plot_frac(table, name, x, y, inv ):

    xlabs = ["Homo","Prim", "Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Vert"]

    if inv == 1: # if you invert axis, make sure everything gets configured
        xlabs = xlabs[::-1] # if you want oldest ages on top.
        table = table.sort_index(ascending =False).sort_index(ascending  =False, axis = 1)


        mask = np.tril(np.ones_like(table, dtype=np.bool))
        outf = f"{RE}{name}_core_der_age_fraction_heatmap_inv.pdf"

    else:
        outf = f"{RE}{name}_core_der_age_fraction_heatmap.pdf"
        mask = np.triu(np.ones_like(table, dtype=np.bool))

    sns.set(font_scale=1.4)
    fig, ax = plt.subplots(figsize=(9,8))

    sns.heatmap(table,
    #mask =mask,
    annot = True,
    linewidths = 1,
    linecolor = "k",
    cmap = "Greys",
    cbar_kws={'label': '% core-derived age pair\nper core age'},
    ax = ax
    )

    ax.set_xticklabels(xlabs)
    ax.set_yticklabels(xlabs)
    ax.set(title = f"{name} % of core-derived age pair",
    ylabel = y,
    xlabel = x)

    plt.savefig(outf, bbox_inches = "tight")
    plt.show()


def pivot_table(x, y, df, val):
    table = pd.pivot(df, index = x, columns = y, values = val)
    table = add_zeros_to_table(table) # add zeros
    table = table.iloc[::-1] # invert

    return table

def plot_OR(OR_log2_table, OR_log10_annot, OR_table, inv):
    xlabs = ["Homo","Prim", "Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Vert"]
    if inv == 1: # if you invert axis, make sure everything gets configured
        xlabs = xlabs[::-1] # if you want oldest ages on top.

        OR_log2_table = OR_log2_table.sort_index(ascending =False).sort_index(ascending  =False, axis = 1)
        OR_log10_annot = OR_log10_annot.sort_index(ascending =False).sort_index(ascending  =False, axis = 1)
        OR_table = OR_table.sort_index(ascending =False).sort_index(ascending  =False, axis = 1)

        mask = np.tril(np.ones_like(OR_table, dtype=np.bool))
        outf = f"{RE}OR_core_der_age_heatmap_annot-log2_inv.pdf"

    else:
        outf = f"{RE}OR_core_der_age_heatmap_annot-log2.pdf"
        mask = np.triu(np.ones_like(OR_log2_table, dtype=np.bool))

    fig, ax = plt.subplots(figsize=(9,8))
    sns.set(font_scale=1.4)

    # plot log2 OR (color) and annotate with log10p asterisks
    sns.heatmap(OR_log2_table,
    cmap = "bwr",
    vmin=-2, vmax=2,
    center = 0,
    linewidths = 1,
    linecolor = 'k',
    annot = False,
    annot_kws={ "color": "k", "ha": 'left',"va": 'top'},
    fmt = '.1d',
    cbar_kws={'label': 'OR\n(log2-scaled)'},


    )

    # plot log2 annotation
    sns.heatmap(OR_table,
    mask = OR_log10_annot == False,
    annot = True,
    annot_kws={"color": "k", "ha": 'center',"va": 'bottom'},
    alpha = 0,
    cbar = False,
    ax = ax
    )


    ax.set(
    title = "FANTOM core-derived age pair",
    xlabel = "core age",
    ylabel = "derived age"
    )

    ax.set_xticklabels(xlabs)
    ax.set_yticklabels(xlabs)

    plt.savefig(outf, bbox_inches = "tight")


#%%


enh = format_syndf(FANTOM)
enh["id"] = "FANTOM"


#%%
# calculate fractions of der/cor combos per core age
corr = get_core_der_frac(enh)

# pivot corr into datatable w/ fraction overlap as value
x,y =  "der_mrca","core_mrca",
corr_table = corr.pivot(index = x, columns = y, values = "core_der_overlap")
corr_table = add_zeros_to_table(corr_table) # add human column to make matrix square

# plot the fraction of core derived pairs taht overlap
name = "FANTOM"
yplot = corr_table.index.name
xplot = corr_table.columns.name
corr_table.iloc[::-1]

plot_frac(corr_table.iloc[::-1], name, xplot, yplot, 1)

#%%


shuf = format_syndf(SHUFF)
shuf["id"] = "SHUFFLE"

shuf_corr = get_core_der_frac(shuf)
shuf_corr_table = shuf_corr.pivot(index = x, columns = y, values = "core_der_overlap")
shuf_corr_table = add_zeros_to_table(shuf_corr_table)

# plot the fraction of core derived pairs taht overlap
name = "shuffle"

yplot = shuf_corr_table.index.name
xplot = shuf_corr_table.columns.name
shuf_corr_table.iloc[::-1]

plot_frac(shuf_corr_table.iloc[::-1], name, xplot, yplot, 1)

#%%

# calculate the odds ratio of observing
fdr_OR = calc_OR(corr, shuf_corr)

fdr_OR.head()
fdr_OR.loc[fdr_OR.comparison_name == "core-0.957_der-0.175"]
fdr_OR.loc[fdr_OR.comparison_name == "core-0.656_der-0.0"]
#%%
# pivot table for OR
x,y =  "der_mrca","core_mrca",

OR_table = pivot_table(x, y, fdr_OR,"OR")
OR_table = OR_table.round(1)
OR_table
# pivot table for log2OR (color )
OR_log2_table = pivot_table(x, y, fdr_OR, "log2")


# pivot table for log10p
OR_log10_annot = pivot_table(x, y, fdr_OR, "reject_null")
OR_log10_annot.replace([-0, 0], False, inplace=True) # replace all the other values with no value

a_table = pivot_table(x, y, fdr_OR,"a")


x,y =  "der_mrca","core_mrca"

#%% plot OR


plot_OR(OR_log2_table, OR_log10_annot, OR_table, 1)



#%% Count table

fig, ax = plt.subplots(figsize=(9,8))
sns.set(font_scale=1.4)
# plot log2 OR
OR_log2_table = OR_log2_table.sort_index(ascending =False).sort_index(ascending  =False, axis = 1)
sns.heatmap(OR_log2_table,
annot = False,
cmap = "Greys",
cbar_kws={'label': 'OR\n(log2-scaled)'},
ax = ax
)
a_table  = a_table.sort_index(ascending =False).sort_index(ascending  =False, axis = 1)
# plot n annotation
a_table = a_table.fillna(0).astype(int)
sns.heatmap(a_table,
#mask = a_table<1,
annot = True, fmt = 'd',
annot_kws={"color": "k", "ha": 'center',"va": 'top'},
alpha = 0,
cbar = False,
ax = ax
)

ax.set(
title = "core-derived age pairing N",
xlabel = "core age",
ylabel = "derived age"
)
xlabs = xlabs[]
ax.set_xticklabels(xlabs)
ax.set_yticklabels(xlabs, rotation = 90)
outf = f"{RE}OR_core_der_age_heatmap_annot-n.pdf"
plt.savefig(outf, bbox_inches = "tight")


#%%
len(enh.loc[enh["core_remodeling"]==1].enh_id.unique())# 10942

len(shuf.loc[shuf["core_remodeling"]==1].enh_id.unique()) # 1126876

#%%
