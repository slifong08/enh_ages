#!/usr/bin/env python
# coding: utf-8

# In[4]:


# 2019-06-03 - created and run on common variants AF > 0.01

##### updates #####

# 2019-06-10
    # AF calculation was not correct. 0.01<= AF <0.5 is the maf, and 0.5=<AF=<1 is the AF
    # SNPs recalculated as maf.
    # Instead of intersecting only common variants (AF >= 0.01), intersect all variants.
    # maf rounded to 7th decimal


# Sarahfong

# Analyze the genomic shuffle of FANTOM eRNA enhancers for breaks v. actual transcribed enhancers.


import glob
import pandas
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os, sys
from scipy import stats
import seaborn as sns
import subprocess

colors = ["faded green", "greyish",  "amber", "dusty purple", "windows blue",]
palette = sns.xkcd_palette(colors)
sns.palplot(palette)

RE = "/dors/capra_lab/projects/enhancer_ages/1000g/results/"

#Fantom
PATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/architecture_coordinates/"
fs = [
    f"{PATH}all_fantom/no-exon_FANTOM_complexcore_age_breaks.bed",
    f"{PATH}all_fantom/no-exon_FANTOM_complexenh_age_breaks.bed",
    f"{PATH}all_fantom/no-exon_FANTOM_derived_age_breaks.bed",
    f"{PATH}all_fantom/no-exon_FANTOM_simple_age_breaks.bed"]

#1000g
THOU_PATH = "/dors/capra_lab/projects/enhancer_ages/1000g/data/maf/"
common_var_beds = glob.glob(f"{THOU_PATH}trimmed.chr*.phase3_shapeit2_mvncall_integrated_v5a.maf.bed")




#%% subtract exons function


def load_syn_gen_bkgd():

    F = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_taxon.bed"
    syngenbkgd = pandas.read_csv(F, sep='\t')
    syngenbkgd[["mrca", "mrca_2"]] = syngenbkgd[["mrca", "mrca_2"]].round(3)

    return syngenbkgd


def intersect_1000G(file_list, common_var_beds, thou_path):

    F_filter = ["complexenh", "simple"]
    outfile_list = []
    for f in file_list:

        if "shuffle" in f:
            sid = (f.split("/")[-1]).split(".")[0] # different str.split strategy for shuffles.
            done = glob.glob("%sshuffle_FANTOM_*_x_trimmed.bed" % (thou_path))
            len_done = len(done) # Check that you haven't done this already.

        else:
            sid = ((f.split("/")[-1]).split("_")[2])
            done = glob.glob("%sFANTOM_*_x_trimmed.bed" % (thou_path))
            len_done = len(done)# Check that you haven't done this already.



        if sid in F_filter and len_done >4:
            for b in common_var_beds:
                chr_num = (b.split("/")[-1]).split(".")[1]
                outfile = "%s%s_x_%s_trimmed.bed" % (thou_path, chr_num, sid)
                bed_cmd = "bedtools intersect -a %s -b %s  -wa -wb> %s" % (fantom, b, outfile)

                subprocess.call(bed_cmd, shell = True)
                outfile_list.append(outfile)
        else:
            outfile_list = done
    return outfile_list


def get_dataframes(file_list, syngenbkgd):

    files_dict= {}
    for file in file_list:
        df = pandas.read_csv(file, sep = '\t', header = None)


        df.columns = ["chr_enh", "start_enh", "end_enh", "mrca", "enh_id",
        "chr_snp", "pos-1", "pos", "qual", "AF", "gwas_id", "overlap"]

        # formatting
        if "shuffle" in file:
            sid = "_".join((file.split('/')[-1]).split("_")[0:3])
            df["code"] = sid.split("_")[2]
            df["dataset"] = "shuffle"
        else:
            sid = "_".join((file.split('/')[-1]).split("_")[0:2])
            df["code"] = sid.split("_")[1]
            df["dataset"] = "FANTOM"

        print(sid)
        files_dict[sid] = df # add to dictionary


    df = pandas.concat(files_dict.values()) # concat enh dictionaries

    df["mrca"] = df["mrca"].round(3)

    df = pandas.merge(df,syngenbkgd, how = "left", on = "mrca") # taxon names

    df["AF"] = df["AF"].astype(float)

    return df


def get_singletons(df):

    sdf = df.loc[(df["AF"] < 0.0004) & (df["overlap"] == 1)] # get only singletons among all variants that overlap enhancers

    sdf_freq = sdf.groupby(["code", "dataset"])["overlap"].sum().reset_index() # sum the overlaps
    sdf_freq.columns = ["code", "dataset", "singleton_count"]

    total_freq = df.groupby(["code"])["overlap"].sum().reset_index() # sum total overlaps
    total_freq.columns = ["code", "df_count"]

    freq_ = pandas.merge(sdf_freq, total_freq, how = "left", on = "code")

    freq_["freq"] = freq_["singleton_count"].divide(freq_["df_count"])

    return freq_


def singleton_foldchange(freq, shuf_freq):
    shuf_freq.columns = ["code", "shuf_dataset", "shuf_singleton_count", "shuf_df_count", "shuf_freq"]
    fold_change = pandas.merge(freq, shuf_freq, how = "left", on = "code")



    fold_change["foldchange"] = fold_change.freq.divide(fold_change.shuf_freq)

    return fold_change


def plot_singletons(plot, singletonfc):

    plot_order =[ "simple", "complexenh", "complexcore" ,"derived"]

    f, (ax1, ax2) = plt.subplots(ncols = 2, figsize = (12, 6))

    x, y = "code", "freq"
    hue = "dataset"
    data = plot

    sns.barplot(
    x = x,
    y =y,
    data = data,
    hue = hue,
    order = plot_order,
    ax = ax1)

    ax1.set(ylabel = "% of arch dataset",
    title = "% Singleton in archicture dataset")

    ax1.set_xticklabels(ax1.get_xticklabels(), rotation =90)
    ax1.legend(bbox_to_anchor=(1.2, 1.0))

    # plot foldchange
    y = "foldchange",
    data = singletonfc

    sns.pointplot(
        x = x,
        y = y,
        data = data,
        order = plot_order,
        ax = ax2,)
        #join = False)
    ax2.set(ylabel = "fold change arch/bkgd",
    title = "fold change\nFANTOM v. shuffle singleton")

    ax2.set_xticklabels(ax2.get_xticklabels(), rotation =90)

    plt.tight_layout()
    outplot = "%ssingletons_FANTOM_arch.pdf" % RE
    plt.savefig(outplot, bbox_inches = "tight")


def plot_arch_maf(cdf):

    colors = ["amber", "dusty purple", "windows blue"]

    palette = sns.xkcd_palette(colors)

    fig, ax = plt.subplots(figsize =(6,6))

    plot_order = ["simple", "complexcore", "derived"]

    x, y = "code", "AF"
    data = cdf.sort_values(by=(["code"]))
    g = sns.barplot( x = x,
    y = y,
    data = data,
    palette = palette,
    order = plot_order) #hue = "code")

    ax.set_xticklabels(ax.get_xticklabels(), rotation = 90, horizontalalignment = "left")
    ax.set(ylabel="MAF", title= "MAF per architecture", xlabel = "", ylim = (0.045, 0.055))
    outplot = "%s1000G_FANTOM_arch_nosingletons.pdf" % RE
    plt.savefig(outplot, bbox_inches = "tight")

    # do stats
    k, kp = stats.kruskal(cdf.AF.loc[cdf.code.str.contains('simple')],
                                   cdf.AF.loc[cdf.code.str.contains('complexcore')],
                                             cdf.AF.loc[cdf.code.str.contains('derived')])
    m, mp = stats.mannwhitneyu(cdf.AF.loc[cdf.code.str.contains('complexcore')],
                                   cdf.AF.loc[cdf.code.str.contains('derived')])
    stat_df = pd.DataFrame({
    "test": ["krusal-wallis", "MWU"],
    "comparison":["simple v.complexcore v. complexderived", "complexcore v. derived"],
    "stat":[k,m],
    "p-value":[kp, mp]
    })

    out_stat_df = "%s1000G_FANTOM_arch_nosingletons.tsv" % RE

    stat_df.to_csv(out_stat_df, sep = '\t', index = False) # write the stats results


def plot_arch_maf_mrca(cdf):
    colors = ["amber","dusty purple", "windows blue", ]

    plot_order = ["simple", "complexcore", "derived"]
    palette = sns.xkcd_palette(colors)
    fig, ax = plt.subplots(figsize =(16,8))

    g = sns.barplot(
    x = "taxon2",
    y = "AF",
    data = cdf.sort_values(by=(["mrca_2"])),
    hue = "code",
    palette = palette,
    hue_order = plot_order) #hue = "code")

    xlabs = ["homo", "prim", "euar", "bore", 'euth', "ther", "mam",
    "amni", "tetr", "vert"]
    ax.set_xticklabels(xlabs, rotation = 90, horizontalalignment = "left")
    ax.set(ylabel="MAF", title = "MAF per architecture")

    outplot = "%s1000G_FANTOM_arch_nosingletons_mrca.pdf" % RE
    plt.savefig(outplot, bbox_inches = "tight")


def get_metrics(cdf, RE, x):

    #counts
    counts = cdf.groupby("code")[x].count()
    counts["metric"] = "counts"

    #medians
    medians = cdf.groupby("code")[x].median()
    medians["metric"] = "median"

    #means
    means = cdf.groupby("code")[x].mean()
    means["metric"] = "mean"

    # file to write dataframe
    outstat_file = f"{RE}1000g_no_singleton_metric_{x}.tsv"
    # concat dataframes
    outstats = pd.concat([counts, medians, means])

    # write dataframe
    outstats.to_csv(outstat_file, sep = '\t', index = False)


def get_snp_density(cdf):

    test = cdf.groupby(["enh_id", "start_enh", "end_enh", "code"])["overlap"].count().reset_index()

    # in this dataframe, the start and end coordinates correspond to the boundaries of syntenic coordinates
    # there may be >1 derived row for an enhancer (e.g. "chr10:101539519-101539906")
    test["seg_len"] = test.end_enh - test.start_enh

    # sum seg_lens and overlapping SNP counts for enhancers that have more than 1 derived entry.
    snp_den = test.groupby(["enh_id", "code"])["overlap", "seg_len"].sum().reset_index()

    snp_den["snp_den"] = snp_den.overlap.divide(snp_den.seg_len) # calculate the SNP density

    return snp_den


#%% # Fantom arch intersection w/ 1000G

outfiles = intersect_1000G(fs, common_var_beds, THOU_PATH)


#%% # shuffle arch intersection w/ 1000G

SHUFPATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/architecture_coordinates/"
shuff_samples = f"{SHUFPATH}shuffle_FANTOM*.bed"
shuff_samples_list = glob.glob(shuff_samples)
shuff_samples_list
shuff_outfiles = intersect_1000G(shuff_samples_list, common_var_beds, THOU_PATH)
shuff_outfiles
#%% # Analysis

syngenbkgd = load_syn_gen_bkgd()

df = get_dataframes(outfiles, syngenbkgd)
shuf = get_dataframes(shuff_outfiles, syngenbkgd)


#%% # Singletons only
#SINGLETONS: 1 mutation/2504 individuals = 0.000399 mutations / individuals.
# get dataframes for singletons
freq = get_singletons(df)

shuf_freq = get_singletons(shuf)


#%% prepare to plot singletons
plot = pandas.concat([shuf_freq, freq])

# get fold change
singletonfc = singleton_foldchange(freq, shuf_freq)

# plot singletons
plot_singletons(plot, singletonfc)

#%% # Exclude singletons

cdf = df.loc[df["AF"] > 0.0004] # FILTER remove singletons + enhancers w/ no variant overlap

get_metrics(cdf, RE, "AF")
plot_arch_maf(cdf)
plot_arch_maf_mrca(cdf)

# # SNP Density x architecture

# In[39]:

snp_den = get_snp_density(cdf)

print(snp_den.groupby("code").median().round(4))

print(snp_den.groupby("code").mean().round(3))

order = ["simple", "complexenh", "complexcore", "derived"]
fig, (ax1, ax2) = plt.subplots(ncols = 2, figsize =(16,8))
sns.boxplot(x = "code", y = "seg_len", data = snp_den,palette = shufpalette,notch = True,
            showfliers = False, order = order, ax= ax1)
ax1.set_title("Length/ architecture")
ax1.set_xticklabels(ax1.get_xticklabels(), rotation = 90, horizontalalignment = "left")

sns.boxplot(x = "code", y = "snp_den", data = snp_den,palette = shufpalette, notch = True,
            showfliers = False, order = order, ax = ax2)
ax2.set_title("SNP density/ architecture")
ax2.set_ylabel("SNP/bp")
ax2.set_xticklabels(ax2.get_xticklabels(), rotation = 90, horizontalalignment = "left")
plt.savefig("%sFANTOM_arch_snp_density_no_singletons.pdf" % RE, bbox_inches = "tight")


# In[70]:


cdf.head()


# In[75]:


test = cdf.groupby(["enh_id", "start_enh", "end_enh", "code", "mrca_2", "taxon2"])["overlap"].count().reset_index()

# in this dataframe, the start and end coordinates correspond to the boundaries of syntenic coordinates
# there may be >1 derived row for an enhancer (e.g. "chr10:101539519-101539906")
test["seg_len"] = test.end_enh - test.start_enh

# sum seg_lens and overlapping SNP counts for enhancers that have more than 1 derived entry.
snp_den = test.groupby(["enh_id", "mrca_2", "taxon2","code"])["overlap", "seg_len"].sum().reset_index()

snp_den["snp_den"] = snp_den.overlap.divide(snp_den.seg_len) # calculate the SNP density

print(snp_den.groupby("code").median().round(4))

print(snp_den.groupby("code").mean().round(3))

order = ["simple", "complexenh", "complexcore", "derived"]
fig, (ax2) = plt.subplots( figsize =(16,8))

sns.boxplot(x = "taxon2", y = "snp_den", data = snp_den.sort_values(by = "mrca_2"),
            palette = shufpalette, notch = True,
            hue = "code",
            showfliers = False, hue_order = order, ax = ax2)
ax2.set_title("SNP density/ architecture")
ax2.set_ylabel("SNP/bp")
ax2.set_xticklabels(ax2.get_xticklabels(), rotation = 90, horizontalalignment = "left")
plt.savefig("%sFANTOM_arch_snp_density_no_singletons.pdf" % RE, bbox_inches = "tight")


# # SNP density x Age

# In[47]:


cdf.code.loc[(cdf.code.str.contains("simple|complexenh"))].unique()


# In[57]:


# select only the simple/complex enhancers
test = cdf.loc[(cdf.code.str.contains("simple|complexenh"))].groupby(["enh_id", "start_enh", "end_enh", "mrca","code"])["overlap"].count().reset_index()

test["seg_len"] = test.end_enh - test.start_enh
test.head()


# In[58]:


# sum seg_lens and overlapping SNP counts for enhancers that have more than 1 derived entry.
snp_den = test.groupby(["enh_id", "code",  "mrca"])["overlap", "seg_len"].sum().reset_index()
snp_den["snp_den"] = snp_den.overlap.divide(snp_den.seg_len) # calculate the SNP density


# In[59]:


snp_den.groupby("mrca").median().round(4)
snp_den["mrca"]= snp_den["mrca"].round(3)


# In[62]:


snp_den.groupby("mrca").mean().round(4)


# In[61]:


snp_den = pandas.merge(snp_den, syngenbkgd, how = "left", on = "mrca")
snp_den.head()


# In[63]:


snp_den.mrca_2.max()


# In[74]:


fig, ax = plt.subplots(figsize =(8,8))
sns.boxplot(x = "taxon2", y = "snp_den", data = snp_den.sort_values(by = "mrca_2"),hue = "code",
            notch = True, showfliers = False)
ax.set_title("SNP density/ architecture\n younger enhancers have more SNPs than older enhancers")
ax.set_xticklabels(ax.get_xticklabels(), rotation = 90, horizontalalignment = "left")
plt.savefig("%sFANTOM_arch_snp_density_mrca_no_singletons.pdf" % RE, bbox_inches = "tight")


# In[56]:


fig, ax = plt.subplots(figsize =(8,6))
plot = cdf.loc[(cdf.code.str.contains("simple")) | cdf.code.str.contains("complexenh")]
group = plot.groupby(["enh_id","AF", "gwas_id", "code"])["mrca"].max().reset_index()
g = sns.pointplot(x = "mrca", y = "AF", data = group.sort_values(by=(["code", "mrca"]))
                , join=False, showfliers = False, hue = "code")
ax.set_xticklabels(ax.get_xticklabels(), rotation = 300, horizontalalignment = "left")
ax.set_ylabel("MAF")
ax.set_title("simple and complex enhancer allele frequencies by age")


# In[52]:


fig, ax = plt.subplots(figsize =(15,8))

g = sns.lineplot(x = "mrca", y = "AF", data = cdf.sort_values(by=(["code", "mrca"])),
#    showfliers = False,
                 hue = "code", palette = shufpalette, hue_order = order)

#ax.set_xticklabels(ax.get_xticklabels(), rotation = 300, horizontalalignment = "left")
ax.set_ylabel("MAF")
#ax.set_xscale("log")
plt.legend(loc='best', bbox_to_anchor=(0.7, 0.5, 0.5, 0.5))
plt.savefig("%sFANTOM_arch_snp_density_mrca_no_singletons.pdf" % RE, bbox_inches = "tight")


# In[53]:


fig, ax = plt.subplots(figsize =(6,6))
sns.distplot(cdf["AF"].loc[df.code.str.contains("simple")],
             label = "simple enh",
             norm_hist=True,
             hist_kws=dict(cumulative=True),
             kde_kws=dict(cumulative=True), )
sns.distplot(cdf["AF"].loc[df.code.str.contains("complexenh")],
             label = "complex enh",
            norm_hist=True,
             hist_kws=dict(cumulative=True),
             kde_kws=dict(cumulative=True), )

mwu, pval = stats.mannwhitneyu(cdf["AF"].loc[cdf.code.str.contains("simple")],                               cdf["AF"].loc[cdf.code.str.contains("complexenh")])
ax.set_xlabel("MAF, mwu = %s, pval = %s"% (mwu, pval))
ax.set_title("simple v. complex enhancer")
ax.set_xlim(-0.01, 0.51)
plt.legend()
plt.savefig("%s1000G_FANTOM_enh_cd.pdf" % RE, bbox_inches = "tight")


# In[130]:


cdf["AF"].groupby(cdf.code).mean().round(3)


# In[131]:


cdf["AF"].groupby(cdf.code).median().round(4)


# In[141]:


fig, ax = plt.subplots(figsize =(6,6))
sns.distplot(cdf["AF"].loc[cdf.code.str.contains("simple")],
             hist = False,
             label = "simple",
             color = "r",

            hist_kws=dict(cumulative=True),
             kde_kws=dict(cumulative=True))
sns.distplot(cdf["AF"].loc[cdf.code.str.contains("complexcore")],
             hist = False,
             label = "complex_core",
            hist_kws=dict(cumulative=True),
             kde_kws=dict(cumulative=True))
sns.distplot(cdf["AF"].loc[cdf.code.str.contains("derived")],
             hist = False,
             label = "complex_derived",
            hist_kws=dict(cumulative=True),
             kde_kws=dict(cumulative=True))
mwu, pval = stats.kruskal(cdf["AF"].loc[cdf.code.str.contains("simple")],                          cdf["AF"].loc[cdf.code.str.contains("complexcore")],                               cdf["AF"].loc[cdf.code.str.contains("derived")])
ax.set_xlabel("MAF\nKruskal = %s, pval = %s"% (mwu, pval))
ax.set_ylabel("Cumulative frequency of dataset")
ax.set_title("1000G simple v. complex core v. derived\nno singletons")
ax.set_xlim(-0.01, 0.51)
plt.legend()


# # common variants maf > 0.01

# In[54]:


common = df.loc[df["AF"]>=0.01]
print("all", len(df), "| common", len(common), "| percent of all variants that are common:", len(common)/len(df))

fig, ax = plt.subplots(figsize =(6,6))

g = sns.boxplot(x = "code", y = "AF", data = common.sort_values(by="mrca"),
    showfliers = False, palette = shufpalette, order = order, notch = True)
#g = sns.swarmplot(x = "mrca", y = "AF", data = df, dodge = True,  hue = "code")
#g.annotate(stats.pearsonr)
ax.set_xticklabels(ax.get_xticklabels(), rotation = 90, horizontalalignment = "left")
ax.set_ylabel("MAF")
ax.set_title("Common variants (MAF >=0.01)")
plt.savefig("%sFANTOM_arch_maf_commonvar.pdf" % RE, bbox_inches = "tight")


# In[142]:


common["AF"].groupby(common.code).mean().round(3)


# In[143]:


common["AF"].groupby(common.code).median().round(3)


# In[55]:


fig, ax = plt.subplots(figsize =(6,6))
sns.distplot(common["AF"].loc[common.code.str.contains("simple")],
             label = "simple enh",
             norm_hist=True,
            hist_kws=dict(cumulative=True),
             kde_kws=dict(cumulative=True))
sns.distplot(common["AF"].loc[common.code.str.contains("complexenh")],
             label = "complex enh",
            norm_hist=True,
            hist_kws=dict(cumulative=True),
             kde_kws=dict(cumulative=True))

mwu, pval = stats.mannwhitneyu(common["AF"].loc[common.code.str.contains("simple")],                               common["AF"].loc[common.code.str.contains("complexenh")])
ax.set_xlabel("Common MAF, mwu = %s, pval = %s"% (mwu, pval))
ax.set_title("simple v. complex enhancer\n Common variants (MAF>=0.01)")
ax.set_xlim(-0.01, 0.51)
plt.legend()
plt.savefig("%s1000G_FANTOM_enh_cd.pdf" % RE, bbox_inches = "tight")


# In[56]:


fig, ax = plt.subplots(figsize =(10,10))
sns.distplot(common["AF"].loc[common.code.str.contains("simple")],
             hist = False,
             label = "simple",
             color = "r",
            hist_kws=dict(cumulative=True),
             kde_kws=dict(cumulative=True))
sns.distplot(common["AF"].loc[common.code.str.contains("complexcore")],
             hist = False,
             label = "complex_core",
            hist_kws=dict(cumulative=True),
             kde_kws=dict(cumulative=True))
sns.distplot(common["AF"].loc[common.code.str.contains("derived")],
             hist = False,
             label = "complex_derived",
            hist_kws=dict(cumulative=True),
             kde_kws=dict(cumulative=True))
mwu, pval = stats.kruskal(common["AF"].loc[common.code.str.contains("simple")],                          common["AF"].loc[common.code.str.contains("complexcore")],                               common["AF"].loc[common.code.str.contains("derived")])
ax.set_xlabel("Kruskal = %s, pval = %s"% (mwu, pval))
ax.set_title("simple v. complex core v. derived\nCommon variants (MAF>=0.01)")
ax.set_xlim(-0.01, 0.51)
plt.legend()
plt.savefig("%s1000G_FANTOM_enh_arch_cd.pdf" % RE, bbox_inches = "tight")


# # Rare variants maf < 0.01
#

# In[57]:


rare = df.loc[(df["AF"]<0.01) &(df["AF"]>0.000399)]
print("all", len(df), "| rare", len(rare), "| percent of all variants that are rare:", len(common)/len(df))

fig, ax = plt.subplots(figsize =(6,6))

g = sns.boxplot(x = "code", y = "AF", data = rare.sort_values(by="mrca"),
    showfliers = False, order = order, palette = shufpalette, notch = True)

ax.set_xticklabels(ax.get_xticklabels(), rotation = 300, horizontalalignment = "left")
ax.set_ylabel("MAF")

ax.set_title("Rare variants (0.000399< MAF <=0.01)")
plt.savefig("%sFANTOM_arch_maf_rarevar.pdf" % RE, bbox_inches = "tight")


# In[150]:


rare["AF"].groupby(rare.code).mean()


# In[151]:


rare["AF"].groupby(rare.code).median()


# In[153]:


fig, ax = plt.subplots(figsize =(6,6))
sns.distplot(rare["AF"].loc[rare.code.str.contains("simple")],
             kde = False,
             bins = 100,
             label = "simple enh",
             norm_hist=True)
sns.distplot(rare["AF"].loc[rare.code.str.contains("complexenh")],

             kde = False,
             bins = 100,
             label = "complex enh",
            norm_hist=True,)

mwu, pval = stats.mannwhitneyu(rare["AF"].loc[rare.code.str.contains("simple")],                               rare["AF"].loc[rare.code.str.contains("complexenh")])
ax.set_xlabel("rare MAF\n mwu = %s, pval = %s"% (mwu, pval))
ax.set_title("simple v. complex enhancer")
#ax.set_xlim(-0.01, 0.51)
plt.legend()


# In[154]:


fig, ax = plt.subplots(figsize =(10,10))
sns.distplot(rare["AF"].loc[rare.code.str.contains("simple")],
             hist = False,
             label = "simple",
             color = "r",
            hist_kws=dict(cumulative=True),
             kde_kws=dict(cumulative=True))
sns.distplot(rare["AF"].loc[rare.code.str.contains("complexcore")],
             hist = False,
             label = "complex_core",
            hist_kws=dict(cumulative=True),
             kde_kws=dict(cumulative=True))
sns.distplot(rare["AF"].loc[rare.code.str.contains("derived")],
             hist = False,
             label = "complex_derived",
            hist_kws=dict(cumulative=True),
             kde_kws=dict(cumulative=True))
mwu, pval = stats.kruskal(rare["AF"].loc[rare.code.str.contains("simple")],                          rare["AF"].loc[rare.code.str.contains("complexcore")],                               rare["AF"].loc[rare.code.str.contains("derived")])
ax.set_xlabel("rare MAF\nKruskal = %s, pval = %s"% (mwu, pval))
ax.set_title("simple v. complex core v. derived")
#ax.set_xlim(-0.01, 0.51)
plt.legend()


# In[ ]:
