#!/usr/bin/env python
# coding: utf-8


import glob
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os, sys
from scipy import stats
import seaborn as sns
import subprocess


BUILD = "hg19"

ENHPATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/all_fantom_enh/ages/"
ENHFILE = "All_FANTOM_x_ENCODE_hg19.bed"
ENHF = os.path.join(ENHPATH, ENHFILE)


PROTBASE = "/dors/capra_lab/data/gene_age/ProteinHistorian_databases/"

PROTDICT = {
"wagner":f"{PROTBASE}PPODv4_PTHR7-OrthoMCL_wagner1.0/HUMAN_PPODv4_PTHR7-OrthoMCL_wagner1.0_age-depth.protein_list",
"dollo":f"{PROTBASE}PPODv4_PTHR7-OrthoMCL_dollo/HUMAN_PPODv4_PTHR7-OrthoMCL_dollo_age-depth.protein_list",
}

RE = "/dors/capra_lab/projects/enhancer_ages/landscape/results/fantom/protein_ages/"

if os.path.exists(RE) == False:
    os.mkdir(RE)


#%% Color palettes


colors = ["faded green", "greyish",  "amber", "dusty purple", "windows blue",]
palette = sns.xkcd_palette(colors)

arch_colors = [ "amber", "dusty purple", "windows blue","greyish"]
arch_palette = sns.xkcd_palette(arch_colors)
sns.palplot(arch_palette)

plt.rcParams.update({'font.size': 15})


#%%


def load_syn_gen_bkgd(build):

    F = f"/dors/capra_lab/projects/enhancer_ages/{build}_syn_taxon.bed"
    syngenbkgd = pd.read_csv(F, sep='\t')
    syngenbkgd[["mrca", "mrca_2"]] = syngenbkgd[["mrca", "mrca_2"]].round(3)

    return syngenbkgd


def get_protein_df(method, prot_dict):

    file = prot_dict[method]

    p_names = ["ensembl_id", "ppod_id", "tf", "tf_rank_age"]
    pDF = pd.read_csv(file,
    sep = '\t',
    comment='#',
    names = p_names)

    print(pDF.head())

    return pDF


def get_enh_df(enh_file):

    cols = ["chr_syn", "start_syn", "end_syn",
    "enh_id","chr", "start", "end",
    "seg_index", "core_remodeling", "core",
    "mrca",
    "chr_tf", "start_tf", "end_tf",
    "tf_id", "peak_len", "tf",
    "cl",  "overlap"
    ]

    df = pd.read_csv(enh_file,
    comment='#',
    sep = '\t',
    header = None)

    df.columns = cols # add column names

    # redo the tf column
    df["tf"] = df["tf_id"].apply(lambda x: x.split("_")[0])

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
    df.loc[df.overlap <6, "tfoverlap_bin"] = 0


    # round the mrca
    df["mrca"] = df["mrca"].round(3)
    return df


def merge_enh_tf_ages(df1, df2, on_x):

    merged = pd.merge(df1, df2, how = "left", on = on_x)

    return merged


def compare_ages(arch1, arch2, age, df):

        #print(f"comparing {arch1} v. {arch2} in core_age for age {age}")

        test = df.loc[df.core_age == age]
        test.tf_rank_age = test.tf_rank_age.astype(int)

        tfage1 = test.loc[test.arch == arch1, "tf_rank_age"]

        if arch1 == "simple": # comparing simple v. complex
            tfage2 = test.loc[test.arch != arch1, "tf_rank_age"]
        else:
            tfage2 = test.loc[test.arch == arch2, "tf_rank_age"]

        result_stat, p = stats.mannwhitneyu(tfage1, tfage2)

        newdf = pd.DataFrame({
        "arch" :[arch1, arch2],
        "col_metric":["core_age", "core_age"],
        "age":[age, age],
        "stat": [result_stat, result_stat],
        "P":[p, p],
        "means": [tfage1.mean(), tfage2.mean()],
        "meds": [tfage1.median(), tfage2.median()],
        "n_tfbs": [len(tfage1), len(tfage2)]
        })

        return newdf


#%%

METHOD = "dollo"

enhdf = get_enh_df(ENHF)
print(enhdf.shape)
pDF = get_protein_df(METHOD , PROTDICT)


#%% merge syn age annotations


syn_gen_bkgd = load_syn_gen_bkgd(BUILD) # get summarized ages

on_x = "mrca"
enhdf = merge_enh_tf_ages(enhdf, syn_gen_bkgd[["mrca", "mrca_2"]], on_x)

enhdf = enhdf.drop(["mrca"], axis = 1) # drop the mrca column. don't need it.
enhdf = enhdf.drop_duplicates()
print(enhdf.shape)
#enhdf.head()


#%% merge enh and protein age file

on_x = "tf"
df = merge_enh_tf_ages(enhdf, pDF, on_x)

#%% include only the enhancers that overlap TFs

df = df.loc[(df.tfoverlap_bin>0)]
df.shape #(381320, 31)


#%% how many tfs have ages? how many do not?

enh_tfcount = enhdf.tf.unique().size # 339

tfs_wo_ages = df.loc[df.tf_rank_age.isna(), "tf"].unique().size # 10 tfs do not have ages

tfs_w_ages = df.loc[~df.tf_rank_age.isna(), "tf"].unique().size # 328 tfs do have ages

df = df.loc[~df.tf_rank_age.isna()] # only consider enhancers w/ TF age overlap

df["tf_rank_age"].unique()
df.shape # (375542, 31)

"""
TFs WITHOUT AGES !

df.loc[df.tf_rank_age.isna(), "tf"].unique()

array([
'ZNF316', 'HNRNPLL', 'AGO1', 'CUX1', 'AGO2',
'MEF2B', 'CCAR2','RBM14', 'ZBTB8A', 'SUPT20H',
], dtype=object)
"""
plotting_pDF = pDF.loc[~pDF.tf_rank_age.isna()]
plotting_pDF["tf_rank_age"].unique()
remove = ['ZNF285', 'SARS2', 'ZNF234', 'RAB4B',
       'SLC25A10', 'ZNF559-ZNF177', 'CGB8'] # remove the genes that dont have ages. These are genes that do not overlap enhancers

plotting_pDF = plotting_pDF.loc[~plotting_pDF["tf_rank_age"].isin(remove)]

remove_enh = list(df["ensembl_id"].unique())

plotting_pDF = plotting_pDF.loc[~plotting_pDF["ensembl_id"].isin(remove_enh)]
allother_prot = plotting_pDF['tf_rank_age'].astype(int)
enh_prot = df['tf_rank_age'].astype(int)

s, p = stats.mannwhitneyu(enh_prot, allother_prot)

print(s, p, allother_prot.mean(), enh_prot.mean())



#%% SIDE QUEST find enhnacers w/ binding sites from multiple ages?

"""
mrca_count = df.groupby(["enh_id", "tf", "tf_rank_age"])["mrca"].count().reset_index()
#sns.histplot(mrca_count.mrca)

mrca_count.loc[mrca_count.mrca == mrca_count.mrca.max()]

enh_id	tf	tf_rank_age	mrca
166187	chr1:234746093-234747674	FOS	9	5
356430	chr9:73033763-73035382	STAT3	9	5

enhdf.loc[enhdf.enh_id == "chr1:234746093-234747674"]
enhdf.loc[enhdf.enh_id == "chr9:73033763-73035382"]
"""


#%% jointplot


df[["tf_rank_age", "mrca_2"]] = df[["tf_rank_age", "mrca_2"]].astype(float)
x = "tf_rank_age"
y = "mrca_2"
data = df.sort_values(by = x)
sns.jointplot(x = x, y = y, data = data)

#%% analyze all ages


simple_tf_ages = df.loc[df.arch == "simple", "tf_rank_age"]
complex_tf_ages = df.loc[df.arch != "simple", "tf_rank_age"]
core_tf_ages = df.loc[df.arch == "complex_core", "tf_rank_age"]
der_tf_ages = df.loc[df.arch == "complex_derived", "tf_rank_age"]

stats.mannwhitneyu(simple_tf_ages, complex_tf_ages) #MannwhitneyuResult(statistic=17478127699.0, pvalue=1.7383826975039448e-06)
stats.mannwhitneyu(core_tf_ages, der_tf_ages) #MannwhitneyuResult(statistic=4290694046.0, pvalue=9.360525587240483e-05)


"""
df.groupby("arch").tf_rank_age.mean()
arch             tf_rank_age (mean)
complex_core       10.952606 # slightly older than simple ages
complex_derived    10.986806 # slightly older than core ages
simple             10.940260 # slightly younger than complex ages.

# problem - confounded by differences in enhancer sequence ages

df.groupby("arch").tf_rank_age.median()
arch            tf_rank_age (median)
complex_core       12.0
complex_derived    12.0
simple             12.0
Name: tf_rank_age, dtype: float64

"""
#%%
fig, ax = plt.subplots()
sns.kdeplot(simple_tf_ages, label = 'simple')
sns.kdeplot(complex_tf_ages, label = "complex")
ax.legend()
plt.show()


#%% analyze TF ages of complex core v. complex derived based on core age.

core_ages = df.groupby("enh_id")["mrca_2"].max().reset_index()
core_ages.columns = ["enh_id", "core_age"]
df = pd.merge(df, core_ages, how = "left", on = "enh_id")
complex = df.loc[df.arch != "simple"].copy() # make a complex-only df
complex.shape # (186221, 32)

complex.enh_id.unique().size # 8732 complex enhancers evaluated.

#%%

results_dict = {}

for age in complex.core_age.unique():

    arch1 = "complex_core"
    arch2 = "complex_derived"

    results = compare_ages(arch1, arch2, age, complex)
    results_dict[age] = results


#%%

rdf = pd.concat(results_dict.values()).sort_values(by="age")
rdf

#%%


colors = ["dusty purple", "windows blue",]
palette = sns.xkcd_palette(colors)

x = "core_age"
y = "tf_rank_age"
data = complex.loc[complex.core_age>0].sort_values(by = "core_age")
hue_order =["complex_core", 'complex_derived']
hue = "arch"

fig, ax = plt.subplots()
sns.barplot(x=x, y=y, data =data, hue = hue, palette = palette, hue_order = hue_order)
ax.set(xlabel = "enh ages", ylabel = "tf ages (mean)", ylim = (10, 12))
xlabs = ["prim", "euar", "bore", "euth", "ther", "mam", "amni", "tert", "vert"]
ax.set_xticklabels(xlabs, rotation = 90)

outf = f"{RE}protein_ages_core_v_der_by_coreage_{METHOD}.pdf"
plt.savefig(outf, bbox_inches ="tight")

#%%


simple_results_dict = {}

for age in df.core_age.unique():

    arch1 = "simple"
    arch2 = "complex"

    results = compare_ages(arch1, arch2, age, df)
    simple_results_dict[age] = results
#%%

rdf = pd.concat(simple_results_dict.values()).sort_values(by="age")
rdf
#%%
colors = ["amber", "dusty purple", "windows blue",]
palette = sns.xkcd_palette(colors)

x = "core_age"
y = "tf_rank_age"
data = df.loc[df.core_age>0].sort_values(by = "core_age")
hue_order =["simple", "complex_core", 'complex_derived']
hue = "arch"

fig, ax = plt.subplots()
sns.barplot(x=x, y=y, data =data, hue = hue, palette = palette, hue_order = hue_order)
ax.set(xlabel = "enh ages", ylabel = "tf ages (mean)", ylim = (10, 12))
xlabs = ["prim", "euar", "bore", "euth", "ther", "mam", "amni", "tert", "vert"]
ax.set_xticklabels(xlabs, rotation = 90)

outf = f"{RE}protein_ages_simple_v_complex_core_der_by_coreage_{METHOD}.pdf"
plt.savefig(outf, bbox_inches ="tight")

#%%

colors = ["amber", "faded green", "dusty purple", "windows blue",]
palette = sns.xkcd_palette(colors)

x = "core_age"
y = "tf_rank_age"
data = df.loc[df.core_age>0].sort_values(by = "core_age")
hue_order =["simple", "complex"]
hue = "overallarch"

df.overallarch.unique()
fig, ax = plt.subplots()
sns.barplot(x=x, y=y, data =data, hue = hue, palette = palette, hue_order = hue_order)
ax.set(xlabel = "enh ages", ylabel = "tf ages (mean)", ylim = (10, 12))
xlabs = ["prim", "euar", "bore", "euth", "ther", "mam", "amni", "tert", "vert"]
ax.set_xticklabels(xlabs, rotation = 90)

outf = f"{RE}protein_ages_simple_v_complex_summed_by_coreage{METHOD}.pdf"
plt.savefig(outf, bbox_inches ="tight")


#%%
METHOD = "wagner"

enhdf = get_enh_df(ENHF)
print(enhdf.shape)
pDF = get_protein_df(METHOD , PROTDICT)


#%% merge syn age annotations


syn_gen_bkgd = load_syn_gen_bkgd(BUILD) # get summarized ages

on_x = "mrca"
enhdf = merge_enh_tf_ages(enhdf, syn_gen_bkgd[["mrca", "mrca_2"]], on_x)

enhdf = enhdf.drop(["mrca"], axis = 1) # drop the mrca column. don't need it.
enhdf = enhdf.drop_duplicates()
print(enhdf.shape)
#enhdf.head()


#%% merge enh and protein age file

on_x = "tf"
df = merge_enh_tf_ages(enhdf, pDF, on_x)

#%% include only the enhancers that overlap TFs

df = df.loc[(df.tfoverlap_bin>0)]
df.shape #(381320, 31)


#%% how many tfs have ages? how many do not?

enh_tfcount = enhdf.tf.unique().size # 339

tfs_wo_ages = df.loc[df.tf_rank_age.isna(), "tf"].unique().size # 10 tfs do not have ages

tfs_w_ages = df.loc[~df.tf_rank_age.isna(), "tf"].unique().size # 328 tfs do have ages

df = df.loc[~df.tf_rank_age.isna()] # only consider enhancers w/ TF age overlap

df["tf_rank_age"].unique()
df.shape # (375542, 31)

"""
TFs WITHOUT AGES !

df.loc[df.tf_rank_age.isna(), "tf"].unique()

array([
'ZNF316', 'HNRNPLL', 'AGO1', 'CUX1', 'AGO2',
'MEF2B', 'CCAR2','RBM14', 'ZBTB8A', 'SUPT20H',
], dtype=object)
"""
plotting_pDF = pDF.loc[~pDF.tf_rank_age.isna()]
plotting_pDF["tf_rank_age"].unique()
remove = ['ZNF285', 'SARS2', 'ZNF234', 'RAB4B',
       'SLC25A10', 'ZNF559-ZNF177', 'CGB8'] # remove the genes that dont have ages. These are genes that do not overlap enhancers

plotting_pDF = plotting_pDF.loc[~plotting_pDF["tf_rank_age"].isin(remove)]

remove_enh = list(df["ensembl_id"].unique())

plotting_pDF = plotting_pDF.loc[~plotting_pDF["ensembl_id"].isin(remove_enh)]
allother_prot = plotting_pDF['tf_rank_age'].astype(int)
enh_prot = df['tf_rank_age'].astype(int)


s, p = stats.mannwhitneyu(enh_prot, allother_prot)

print(s, p, allother_prot.mean(), enh_prot.mean())
#3200212123.0 6.859196029120678e-55 8.610186199342825 9.190420778501473




#%% SIDE QUEST find enhnacers w/ binding sites from multiple ages?

"""
mrca_count = df.groupby(["enh_id", "tf", "tf_rank_age"])["mrca"].count().reset_index()
#sns.histplot(mrca_count.mrca)

mrca_count.loc[mrca_count.mrca == mrca_count.mrca.max()]

enh_id	tf	tf_rank_age	mrca
166187	chr1:234746093-234747674	FOS	9	5
356430	chr9:73033763-73035382	STAT3	9	5

enhdf.loc[enhdf.enh_id == "chr1:234746093-234747674"]
enhdf.loc[enhdf.enh_id == "chr9:73033763-73035382"]
"""


#%% jointplot


df[["tf_rank_age", "mrca_2"]] = df[["tf_rank_age", "mrca_2"]].astype(float)
x = "tf_rank_age"
y = "mrca_2"
data = df.sort_values(by = x)
sns.jointplot(x = x, y = y, data = data)

#%% analyze all ages


simple_tf_ages = df.loc[df.arch == "simple", "tf_rank_age"]
complex_tf_ages = df.loc[df.arch != "simple", "tf_rank_age"]
core_tf_ages = df.loc[df.arch == "complex_core", "tf_rank_age"]
der_tf_ages = df.loc[df.arch == "complex_derived", "tf_rank_age"]

stats.mannwhitneyu(simple_tf_ages, complex_tf_ages) #MannwhitneyuResult(statistic=17532673679.0, pvalue=0.0016993555652411764)
stats.mannwhitneyu(core_tf_ages, der_tf_ages) #MannwhitneyuResult(statistic=4296390663.0, pvalue=0.0006703023245730718)


"""
df.groupby("arch").tf_rank_age.mean()
arch             tf_rank_age (mean)

complex_core       9.190155 # slightly younger than derived
complex_derived    9.228003 # slightly older than cores
simple             9.172469 # slighly younger than complex
Name: tf_rank_age, dtype: float64

# problem - confounded by differences in enhancer sequence ages

df.groupby("arch").tf_rank_age.median()
arch            tf_rank_age (median)
complex_core       9.0
complex_derived    9.0
simple             9.0
Name: tf_rank_age, dtype: float64

"""
#%%
fig, ax = plt.subplots()
sns.kdeplot(simple_tf_ages, label = 'simple')
sns.kdeplot(complex_tf_ages, label = "complex")
ax.legend()
plt.show()


#%% analyze TF ages of complex core v. complex derived based on core age.

core_ages = df.groupby("enh_id")["mrca_2"].max().reset_index()
core_ages.columns = ["enh_id", "core_age"]
df = pd.merge(df, core_ages, how = "left", on = "enh_id")
complex = df.loc[df.arch != "simple"].copy() # make a complex-only df
complex.shape # (186221, 32)

complex.enh_id.unique().size # 8732 complex enhancers evaluated.

#%%

results_dict = {}

for age in complex.core_age.unique():

    arch1 = "complex_core"
    arch2 = "complex_derived"

    results = compare_ages(arch1, arch2, age, complex)
    results_dict[age] = results


#%%

rdf = pd.concat(results_dict.values()).sort_values(by="age")
rdf

#%%

set(complex.arch)
#%%
colors = ["dusty purple", "windows blue",]
palette = sns.xkcd_palette(colors)

x = "core_age"
y = "tf_rank_age"
data = complex.loc[complex.core_age>0].sort_values(by = "core_age")
hue_order =["complex_core", 'complex_derived']
hue = "arch"

fig, ax = plt.subplots()
sns.barplot(x=x, y=y, data =data, hue = hue, palette = palette, hue_order = hue_order)
ax.set(xlabel = "enh ages", ylabel = "tf ages (mean)", ylim = (8, 10))
xlabs = ["prim", "euar", "bore", "euth", "ther", "mam", "amni", "tert", "vert"]
ax.set_xticklabels(xlabs, rotation = 90)

outf = f"{RE}protein_ages_core_v_der_by_coreage_{METHOD}.pdf"
plt.savefig(outf, bbox_inches ="tight")
#%%


simple_results_dict = {}

for age in df.core_age.unique():

    arch1 = "simple"
    arch2 = "complex"

    results = compare_ages(arch1, arch2, age, df)
    simple_results_dict[age] = results
#%%

rdf = pd.concat(simple_results_dict.values()).sort_values(by="age")
rdf
#%%
colors = ["amber", "dusty purple", "windows blue",]
palette = sns.xkcd_palette(colors)

x = "core_age"
y = "tf_rank_age"
data = df.loc[df.core_age>0].sort_values(by = "core_age")
hue_order =["simple", "complex_core", 'complex_derived']
hue = "arch"

fig, ax = plt.subplots()
sns.barplot(x=x, y=y, data =data, hue = hue, palette = palette, hue_order = hue_order)
ax.set(xlabel = "enh ages", ylabel = "tf ages (mean)", ylim = (8, 10))
xlabs = ["prim", "euar", "bore", "euth", "ther", "mam", "amni", "tert", "vert"]
ax.set_xticklabels(xlabs, rotation = 90)

outf = f"{RE}protein_ages_simple_v_complex_core_der_by_coreage_{METHOD}.pdf"
plt.savefig(outf, bbox_inches ="tight")
#%%

colors = ["amber", "faded green", "dusty purple", "windows blue",]
palette = sns.xkcd_palette(colors)

x = "core_age"
y = "tf_rank_age"
data = df.loc[df.core_age>0].sort_values(by = "core_age")
hue_order =["simple", "complex"]
hue = "overallarch"

df.overallarch.unique()
fig, ax = plt.subplots()
sns.barplot(x=x, y=y, data =data, hue = hue, palette = palette, hue_order = hue_order)
ax.set(xlabel = "enh ages", ylabel = "tf ages (mean)", ylim = (8, 10))
xlabs = ["prim", "euar", "bore", "euth", "ther", "mam", "amni", "tert", "vert"]
ax.set_xticklabels(xlabs, rotation = 90)

outf = f"{RE}protein_ages_simple_v_complex_summed_by_coreage{METHOD}.pdf"
plt.savefig(outf, bbox_inches ="tight")
