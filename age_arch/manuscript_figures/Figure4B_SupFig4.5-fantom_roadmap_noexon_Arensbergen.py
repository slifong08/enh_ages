#!/usr/bin/env python
# coding: utf-8

# In[2]:


import glob
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator
import numpy as np
import os, sys
import pybedtools as pb
from scipy import stats
import seaborn as sns
import statsmodels.api as sm

# sns colors
arch_colors = [ "amber", "dusty purple", "windows blue","greyish"]
arch_palette = sns.xkcd_palette(arch_colors)

colors = [ "amber", "faded green", "dusty purple", "windows blue","greyish"]
palette = sns.xkcd_palette(colors)
sns.palplot(arch_palette)

# sns graphing preferences
sns.set(color_codes=True)
sns.set(font_scale=1.5)
sns.set_style("white")
sns.despine(bottom=True, left=True)

import datetime
last_run = datetime.datetime.now()
today = (datetime.date.today())
print("last run", datetime.datetime.now())

RE = "/dors/capra_lab/projects/enhancer_ages/arensbergen2019/results/"


# In[74]:


RUN_BED = 1
RUN_OR = 1


# In[ ]:


sid_dict = {"UBERON_0002107": "/dors/capra_lab/projects/enhancer_ages/fantom/data/download/UBERON_0002107_liver_expressed_enhancers/UBERON_0002107_liver_expressed_enhancers_enh_ages_enh_age_arch_summary_matrix.bed",
            "CL_0000094": "/dors/capra_lab/projects/enhancer_ages/fantom/data/download/CL_0000094_granulocyte_expressed_enhancers/CL_0000094_granulocyte_expressed_enhancers_enh_ages_enh_age_arch_summary_matrix.bed",
            "E123": "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/Hsap_H3K27ac_plus_H3K4me3_minus_E123/trimmed/no-exon_trimmed-310-E123.bed",
            "E118": "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/Hsap_H3K27ac_plus_H3K4me3_minus_E118/trimmed/no-exon_trimmed-310-E118.bed",
            "all_fantom" : "/dors/capra_lab/projects/enhancer_ages/fantom/data/trimmed_all_unique_fantom_enh_112_tissue.bed"
           }

arens_path = "/dors/capra_lab/data/mpra/arensbergen18/"
outdata_path = "/dors/capra_lab/projects/enhancer_ages/arensbergen2019/data/"

for sid, f in sid_dict.items():

    print(sid)

    enh = pb.BedTool(f).sort()

    outfile = "%sarens_%s.bed" % (outdata_path, sid)

    if RUN_BED ==1:

        arens = pb.BedTool("%sSuRE_SNP_table_181029_hg19.bed" % (arens_path)).sort()
        arens.intersect(enh, wao = True).saveas(outfile)


# In[43]:


fs = glob.glob("%sarens_*.bed" % (outdata_path))
fs


# In[49]:


#filter for the files you want to analyze?

select_sids =["UBERON_0002107", "E123", "CL_0000094", "E118", "all_fantom"]

fsnew = {}

for f in fs:

    f_sid = "_".join(((f.split("/")[-1]).split(".")[0]).split("_")[1:])
    if f_sid in select_sids:

        fsnew[f_sid] = f
fsnew.items()


# In[130]:


def format_df(sid, df):

    fantom_sids = ["UBERON_0002107", "CL_0000094"]
    roadmap_sids = ["E123", "E118"]
    # name the columns. This is annoying bc I made different file formats.
    if sid in fantom_sids:
        df = df[[0, 1,2,3,4,5,6,7,8,9,10,11,12,13, 17, 19, 20, 21, 22, 23, 32]]

        columns = ["chr","start_snp", "end_snp",
                   "SNP_ID","ref.element.count", "alt.element.count",
                   "k562.ref.mean","k562.alt.mean", "k562.wilcoxp","k562.wilcoxp_random",
                   "hepg2.ref.mean", "hepg2.alt.mean","hepg2.wilcoxp","hepg2.wilcoxp_random",
               "enh_id", "sid", "core_remodeling","arch", "seg_index", "mrca","overlap"]
        df.columns = columns
        df["core_remodeling"] =0

    elif sid in roadmap_sids:
        df = df[[0, 1,2,3,4,5,6,7,8,9,10,11,12,13, 17, 19, 20, 21, 30]]

        columns = ["chr","start_snp", "end_snp", "SNP_ID","ref.element.count",
              "alt.element.count","k562.ref.mean","k562.alt.mean",
              "k562.wilcoxp","k562.wilcoxp_random","hepg2.ref.mean",
              "hepg2.alt.mean","hepg2.wilcoxp","hepg2.wilcoxp_random",
                "enh_id", "arch",  "seg_index","mrca",
               "overlap"]

        df.columns = columns
        df["core_remodeling"] =0
        df.loc[df.seg_index ==1, "core_remodeling"] = 1


    elif sid == "all_fantom":
        columns = ["chr","start_snp", "end_snp",
           "SNP_ID","ref.element.count","alt.element.count",
           "k562.ref.mean","k562.alt.mean", "k562.wilcoxp","k562.wilcoxp_random",
           "hepg2.ref.mean", "hepg2.alt.mean","hepg2.wilcoxp","hepg2.wilcoxp_random",
           "chr_enh", "start_enh", "end_enh",
             "enh_id", "core_remodeling", "mrca", "arch", "dataset", "dataset2", "overlap"]
        df.columns = columns
        df["seg_index"] = 1


    df = df.drop_duplicates() # drop any synteny linked duplicates.

    # replace all the missing data with a 1
    df[["k562.ref.mean","k562.alt.mean",
          "k562.wilcoxp","k562.wilcoxp_random","hepg2.ref.mean",
          "hepg2.alt.mean","hepg2.wilcoxp","hepg2.wilcoxp_random"]] = df[["k562.ref.mean","k562.alt.mean",
          "k562.wilcoxp","k562.wilcoxp_random","hepg2.ref.mean",
          "hepg2.alt.mean","hepg2.wilcoxp","hepg2.wilcoxp_random"]].replace(".",1)

    # float columns
    df[["k562.ref.mean","k562.alt.mean", "k562.wilcoxp","hepg2.ref.mean",
                  "hepg2.alt.mean", "hepg2.wilcoxp"]]\
    = df[["k562.ref.mean","k562.alt.mean", "k562.wilcoxp","hepg2.ref.mean",
                  "hepg2.alt.mean", "hepg2.wilcoxp"]].astype(float)

    # reduce dataframe
    df2 = df[["enh_id",
               "k562.ref.mean", "k562.alt.mean","k562.wilcoxp", "core_remodeling",
               "hepg2.ref.mean", "hepg2.alt.mean", "hepg2.wilcoxp", "overlap", "seg_index"]]\
    .drop_duplicates()

    # measure delta activity
    df2["delta_k562"] = (df2["k562.ref.mean"] - df2["k562.alt.mean"])
    df2["delta_hepg2"] = (df2["hepg2.ref.mean"] - df2["hepg2.alt.mean"])

    # measure normalized delta activity to reference mean
    df2["norm_delta_k562"] = (df2["delta_k562"]).divide(df2["k562.ref.mean"])
    df2["norm_delta_hepg2"] = (df2["delta_hepg2"]).divide(df2["hepg2.ref.mean"])

    # replace infinite values with 0
    df2 = df2.replace([np.inf, -np.inf], 0)

    # created significant activity 0|1 binaries
    # based on nominal p-value = 0.05 or 5% FDR p-value cut off determined in arensbergen 2019 publication
    df2["sig_k562_05"],  df2["sig_hepg2_05"], df2["sig_k562_fdr05"], df2["sig_hepg2_fdr05"] = 0, 0, 0, 0

    df2.loc[df2["k562.wilcoxp"]<0.05, "sig_k562_05"] =1
    df2.loc[df2["hepg2.wilcoxp"]<0.05, "sig_hepg2_05"] =1

    df2.loc[df2["k562.wilcoxp"]<0.00173121, "sig_k562_fdr05"] =1
    df2.loc[df2["hepg2.wilcoxp"]<0.006192715, "sig_hepg2_fdr05"] =1

    return df2


# In[69]:


def bootstrapCI(fold_changes_list): # bootstrap C.I.s

    n = len(fold_changes_list)

    xbar = np.mean(fold_changes_list) # get the real median fold-change from 500 shuffles

    nboot = 10000 # resample 10000 times
    val = 0
    bs_means = []

    while val < nboot:

        bs_dist = np.random.choice(fold_changes_list, replace = True, size = n)
        bsmean = np.mean(bs_dist)
        bs_means.append(bsmean)
        val +=1

    bs = pandas.DataFrame(data = bs_means, index = np.arange(nboot), columns = ["bs_means"])

    bs["deltas"] = bs.bs_means - xbar

    bs = bs.sort_values(by = "deltas", ascending= False)

    low = bs.deltas.quantile(0.025)
    high = bs.deltas.quantile(0.975)
    print(low, high)

    ci = xbar - [high, low]
    return ci


# In[159]:


def get_OR(obs, comparison, cell_model):
    print(obs, comparison, cell_model)
    odds, p = stats.fisher_exact(obs)
    table = sm.stats.Table2x2(obs) # get confidence interval
    odds_ci = table.oddsratio_confint()


    results = pd.DataFrame({"sid": [sid],
            "cell_model": [cell_model],
            "nominal_p":[nom_p],
            "a": [obs[0][0]],
            "b": [obs[0][1]],
            "c": [obs[1][0]],
            "d": [obs[1][1]],
            "OR": [odds],
            "p": [p],
            "test":[comparison],
            "ci_lower" :[odds_ci[0]],
            "ci_upper" :[odds_ci[1]]})

    return results


# In[170]:


models = ["K562", "HEPG2"]
nom_p = "FDR05"

ci = {}
all_results = {} # collect all results for both K562 and HepG2 cell lines

all_results.keys()


# In[171]:


if RUN_OR ==1:

    for sid, f in fsnew.items():
        print(sid, f)

        f_results = {}

        key = str(sid+"FDR05")


        if key not in all_results.keys():

            print(sid)

            f_results = {} # collect results per enhancer dataset


            df = pd.read_csv(f, sep = '\t', header = None) # open the file

            df2 = format_df(sid, df) # format dataframe

            # reassign enhancer architectures (aka core_remodeling). If not enhancer, arch = -1
            if sid != "all_fantom":
                df2.loc[df2.seg_index !="1", "core_remodeling"] = 1
                df2.loc[df2.seg_index =="1", "core_remodeling"] = 0
                df2.loc[df2.seg_index ==".", "core_remodeling"] = -1
            print(df2.core_remodeling.unique())

            for CELL_MODEL in models:

                print(CELL_MODEL)

                if CELL_MODEL == 'K562': # create value to pick cell-model matching significance column
                    sig = "sig_k562_fdr05"

                elif CELL_MODEL == 'HEPG2':
                    sig = "sig_hepg2_fdr05"

                ### Get numbers for 2x2 table ###

                # Simple
                simple_sig = df2.loc[(df2.core_remodeling ==0) &
                                              (df2.overlap == 1) &
                                              (df2[sig]==1)]["enh_id"].count()
                simple_total = df2.loc[(df2.core_remodeling ==0) & df2.overlap ==1]["enh_id"].count()
                simple_nonsig = df2.loc[(df2.core_remodeling ==0) &
                                  (df2.overlap == 1) &
                                  (df2[sig]==0)]["enh_id"].count()
                # Complex
                complex_sig = df2.loc[(df2.core_remodeling ==1) &
                                      (df2.overlap == 1) &
                                      (df2[sig]==1)]["enh_id"].count()
                complex_total = df2.loc[(df2.core_remodeling ==1) & df2.overlap ==1]["enh_id"].count()
                complex_nonsig = df2.loc[(df2.core_remodeling ==1) &
                                  (df2.overlap == 1) &
                                  (df2[sig]==0)]["enh_id"].count()

                # bkgd (w/o enhancer overlap)
                bkdg_sig = len(df2.loc[(df2.core_remodeling ==-1) &
                                     (df2[sig]==1) &
                                     (df2.overlap == 0)]) # all sig. results that don't overlap fantom enhancers
                bkdg_nonsig = len(df2.loc[(df2.core_remodeling ==-1)&
                                     (df2[sig]==0) &
                                     (df2.overlap == 0)])# all nonsig. results that don't overlap fantom enhancers


                # Test 1 - sig simple v. all background
                #[[simple sig, simple not sig], [allbkgd sig(minus enh), allbkgd not sig(minus enh)]]


                COMPARISON1 = "sig simple v. bkgd"
                obs1 = np.array([[simple_sig, simple_nonsig], [bkdg_sig, bkdg_nonsig]]) # make a contigency table

                simple_v_bkgd = get_OR(obs1, COMPARISON1, CELL_MODEL) # calculate OR

                key = CELL_MODEL + str(nom_p) + COMPARISON1
                f_results[key] = simple_v_bkgd


                # Test 2 - sig complex v. all background
                #[[complex sig, complex not sig], [allbkgd sig(minus enh), allbkgd not sig(minus enh)]]


                COMPARISON2 = "sig complex v. bkgd"
                obs2 = np.array([[complex_sig, complex_nonsig], [bkdg_sig, bkdg_nonsig]]) # make a contigency table

                complex_v_bkgd = get_OR(obs2, COMPARISON2, CELL_MODEL)# calculate OR

                key = CELL_MODEL + str(nom_p) + COMPARISON2
                f_results[key] = complex_v_bkgd


                # Test 3 - sig simple v. Signsigificant complex
                #[[complex sig, complex not sig], [allbkgd sig(minus enh), allbkgd not sig(minus enh)]]


                COMPARISON3 = "sig simple v. sig complex"
                obs3 = np.array([[simple_sig, simple_nonsig], [complex_sig, complex_nonsig]]) # make a contigency table

                simple_v_complex = get_OR(obs3, COMPARISON3, CELL_MODEL)# calculate OR

                key = CELL_MODEL + str(nom_p) + COMPARISON3
                print(key)
                f_results[key] = simple_v_complex


            concat_results = pd.concat(f_results.values()) # concat all results for this cell line
            key = sid + str(nom_p)
            all_results[key]= concat_results

    df = pd.concat(all_results.values())
    df.sort_values(by ="sid")
    df.to_csv("%sall_trimmed_or.tsv"% (outdata_path), sep = '\t', header= True, index = False)


# In[191]:


df.sid.unique()


# In[196]:


if RUN_OR ==0:
    df = pd.read_csv("%sall_trimmed_or.tsv"% (outdata_path), sep = '\t')

# match cell-specific enhancers with cell model background tested in MPRA

df["cell_line_match"] = 0 # create dummy column

hepg2_sids = ["UBERON_0002107", "E118", "all_fantom"] # hepg2/liver facet ids

k562_sids = ["CL_0000094", "E123", "all_fantom"] # k562/granulocyte facet ids

df["cell_line_match"].loc[(df.cell_model =="HEPG2") & (df.sid.isin(hepg2_sids))] = 1
df["cell_line_match"].loc[(df.cell_model =="K562") & (df.sid.isin(k562_sids))] = 1

df["-log10p"] = -1*np.log10(df.p) # -log10p

df["log2or"] = np.log2(df.OR) # log2(OR)

# dataset names
df['dataset'] = "Fantom_all"
df.loc[df.sid.str.contains("E1"), "dataset"] = "Roadmap-Matched Cell line"
df.loc[df.sid.str.contains("000"), "dataset"] = "Fantom-Matched Tissue"
df["test2"] = df["cell_model"] + "-" + df["test"] # label for plotting Figure 4a


# In[200]:


data = df.loc[(df.cell_line_match==1) & (df.nominal_p.astype(str) != "0.05")]
data.sort_values(by = "p")


# In[201]:


# No background overlap is evaluated here. We don't have background expectation for simple v. complex allele enrichment

plot = data.loc[data.sid != "all_fantom"] # focus on allele overlap between simple and complex only.

plot["SUMab"] = plot.a + plot.b # total alleles in simple enhancers

plot["SUMcd"] = plot.c + plot.d # total alleles in complex enhancers
plot["total"] = plot.a + plot.b + plot.c + plot.d # total alleles in either enhancer

plot["FRACa"] = plot.a.divide(plot.SUMab) # fraction of sig alleles in simple enhancers
plot["FRACab"] = plot.SUMab.divide(plot.total) # fraction of simple alleles in total allele dataset
plot["FRACc"] = plot.c.divide(plot.SUMcd) # fraction of sig alleles in complex enhancers
plot["FRACcd"] = plot.SUMcd.divide(plot.total) # fraction of complex alleles in total allele dataset


# In[202]:


plot.loc[(plot.cell_model == "K562") & (plot.test.str.contains("bkgd"))]


# In[205]:


plot.loc[(plot.sid ==sid) & (plot.test.str.contains("bkgd")), ["p", "OR"]]


# In[207]:


hue_order = [ "CL_0000094", "UBERON_0002107","E123","E118", ]
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize = (12,12))
val = 0
ax_list = [ax1, ax2, ax3, ax4]
for sid in hue_order:
    ax = ax_list[val]
    sns.barplot(x = "cell_model", y = "FRACa",
            data = plot.loc[(plot.sid ==sid) & (plot.test.str.contains("bkgd"))], hue = "test",
             ax = ax,palette = palette)
    bkgd =  round(plot.loc[(plot.sid ==sid) & (plot.test.str.contains("bkgd")), "FRACc"].unique()[0], 3)
    ax.axhline(plot.loc[(plot.sid ==sid) & (plot.test.str.contains("bkgd")), "FRACc"].unique()[0], color = "k", ls = "--")
    ax.set_title(sid)
    ax.legend().remove()

    OR = "OR =" + plot.loc[(plot.sid ==sid) & (plot.test.str.contains("sig simple v. sig complex")), "OR"].round(3).astype(str).item()
    pval = "p =" + plot.loc[(plot.sid ==sid) & (plot.test.str.contains("sig simple v. sig complex")), "p"].round(3).astype(str).item()
    ax.set(xlabel = (OR, pval, bkgd), ylabel= "Frac of sig. MPRA alleles")
    val+=1
plt.tight_layout()
plt.savefig("%sFigS4.5_Arens_pcutoff-arens_Simple_frac_no_exon_trimmed.pdf" %RE)


# In[210]:


data.sid.unique()


# In[211]:


# No background overlap is evaluated here. We don't have background expectation for simple v. complex allele enrichment

plot_fantom = data.loc[data.sid == "all_fantom"] # focus on allele overlap between simple and complex only.

plot_fantom["SUMab"] = plot_fantom.a + plot_fantom.b # total alleles in simple enhancers

plot_fantom["SUMcd"] = plot_fantom.c + plot_fantom.d # total alleles in complex enhancers
plot_fantom["total"] = plot_fantom.a + plot_fantom.b + plot_fantom.c + plot_fantom.d # total alleles in either enhancer

plot_fantom["FRACa"] = plot_fantom.a.divide(plot_fantom.SUMab) # fraction of sig alleles in simple enhancers
plot_fantom["FRACab"] = plot_fantom.SUMab.divide(plot_fantom.total) # fraction of simple alleles in total allele dataset
plot_fantom["FRACc"] = plot_fantom.c.divide(plot_fantom.SUMcd) # fraction of sig alleles in complex enhancers
plot_fantom["FRACcd"] = plot_fantom.SUMcd.divide(plot_fantom.total) # fraction of complex alleles in total allele dataset


# In[212]:


plot_fantom.loc[(plot_fantom.test.str.contains("bkgd")), ["p", "OR", "test", "cell_model"]]


# In[213]:


plot_fantom.loc[(plot_fantom.test.str.contains("bkgd"))]


# In[214]:


colors = [ "amber", "faded green", "greyish", "amber", "faded green", "greyish"]
palette = sns.xkcd_palette(colors)
sns.palplot(palette)


# In[216]:


hue_order = ['all_fantom']
fig, (ax1, ax2,) = plt.subplots(ncols = 2, figsize = (12,6))
val = 0
ax_list = [ax1, ax2,]
for sid in hue_order:
    for cell_model in plot_fantom.cell_model.unique():
        test = plot_fantom.loc[(plot_fantom.sid ==sid) & (plot_fantom.cell_model == cell_model)]
        ax = ax_list[val]
        sns.barplot(x = "cell_model", y = "FRACa",
                data =test.loc[test.test.str.contains("bkgd")],
                    hue = "test",
                 ax = ax,palette = palette)
        bkgd =  round( test["FRACc"].unique()[0], 3)
        ax.axhline(test["FRACc"].unique()[0], color = "k", ls = "--")
        ax.set_title(sid)
        ax.legend().remove()

        OR = "OR =" + test.loc[(test.test.str.contains("sig simple v. sig complex")), "OR"].round(3).astype(str).item()
        pval = "p =" + test.loc[ (test.test.str.contains("sig simple v. sig complex")), "p"].round(3).astype(str).item()
        ax.set_xlabel( (OR, pval, bkgd))
        ax.set_ylabel("Frac of sig. MPRA alleles")
        ax.set_ylim(0,0.05)
        val+=1
plt.tight_layout()
plt.savefig("%sFig4b_Arens_trimmed_fantom_enhancers.pdf" %RE)


# In[218]:


fig, (ax, ax2) = plt.subplots(ncols = 2, figsize = (12,6))



sns.barplot(x = "test2", y = "FRACa",
            data = df.loc[(df.cell_model == "K562") & (df.test.str.contains("bkgd"))],
            ax = ax, palette = palette)

ax.axhline(df.loc[(df.cell_model == "K562") & (df.test.str.contains("bkgd")), "FRACc"].unique()[0], color = "k")
ax.set_ylim(0,0.05)
sns.barplot(x = "test2", y = "FRACa",
            data = df.loc[(df.cell_model == "HEPG2")& (df.test.str.contains("bkgd"))],
            ax = ax2, palette = palette)

ax2.axhline(df.loc[(df.cell_model == "HEPG2") & (df.test.str.contains("bkgd")), "FRACc"].unique()[0], color = "k") # dataset background
ax2.set_ylim(0,0.05)
plt.tight_layout()
plt.savefig("%sfig4b-trimmed_Arens_pcutoff-arens_FRAC_barplot_nomp.pdf" %RE, bbox_inches = "tight")


# In[ ]:


df.loc[(df.cell_model == "K562") & (df.test.str.contains("bkgd")), "FRACc"].unique()[0]


# In[ ]:


df


# In[ ]:
