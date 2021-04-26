#!/usr/bin/env python
# coding: utf-8

# In[1]:


import datetime
import glob
import matplotlib.pyplot as plt
import numpy as np
import os, sys
import pandas as pd
import pybedtools as pb
import seaborn as sns
from scipy import stats
from statsmodels.stats import multitest
import statsmodels.api as sm
import subprocess

today = datetime.date.today()
RE = "/dors/capra_lab/projects/enhancer_ages/fantom/results/te/310-"

print(today)
get_ipython().run_line_magic('matplotlib', 'inline')

colors = ["amber", "dusty purple","windows blue"]
palette = sns.xkcd_palette(colors)
sns.palplot(palette)

es_colors = ["amber", "faded green"]
es_pal = sns.xkcd_palette(es_colors)
sns.palplot(es_pal)
plt.rcParams.update({'font.size': 15})


# In[2]:


TRIMMED = 1 # run syntenic block file or enhancer file?
BED_INTERSECTION = 0

te_path = "/dors/capra_lab/users/abraha1/projects/transposable_elements/data/"
tefam = "%shg19_TE_counts_wlin.txt" % te_path
dfam_file="%sdfam/hg38_dfam.nrph.hits_genomicCoordinates_tidy.bed" % te_path
repeatmasker = "%sfiltered_formated_hg19fromhg38-TE_coords.tsv" % te_path

syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2", "mya", "mya2"]] # whittle down the df
syn_gen_bkgd.head()


# In[3]:


# Fantom x TE intersection

fantom_path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/"
outpath= "/dors/capra_lab/projects/enhancer_ages/fantom/data/te/"

if TRIMMED ==1:
    fantom_file = "%strimmed-310_all_unique_fantom_erna_112_tissue.bed" % fantom_path
else:
    fantom_file = "%sall_unique_fantom_erna_112_tissue.bed" % fantom_path

sid = (fantom_file.split("/")[-1]).split(".")[0]
print(sid)

outf = "%s%s_x_te.bed" % (outpath, sid)

# bedtools intersect
if BED_INTERSECTION == 1:

    f = pb.BedTool(fantom_file).sort()
    r = pb.BedTool(repeatmasker).sort()

    f.intersect(r, wao = True).saveas(outf)


# In[4]:


# shuffle x TE intersection

BED_INTERSECTION = 0
shuf_path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/shuffle/te/"

if TRIMMED == 1:
#    shuf_path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/shuffle/breaks/"

    shuffles = glob.glob("%strimmed-310_shuf-*_age_breaks.bed"% shuf_path)
else:
    shuffle_path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/shuffle/breaks/"

    shuffles = glob.glob("%sshuf-*.bed"% shuffle_path)

for shuf in shuffles:
    shuf_id = (shuf.split("/")[-1]).split(".")[0]
    if "2501x" not in shuf_id:

        # bedtools intersect
        if BED_INTERSECTION == 1:
            s = pb.BedTool(shuf).sort()
            r = pb.BedTool(repeatmasker).sort()
            shuf_outf = "%s%s_x_te.bed" % (shuf_path, shuf_id)

            s.intersect(r, wao = True).saveas(shuf_outf)


# In[5]:


def format_df(df, data_id, trimmed):
    if trimmed==1:
        df = df[[0,1,2,3,4,5,6,7,8,9,10,11]]
        df.columns = ["chr_enh", "start_enh", "end_enh", "old_len", "enh_id", "mrca","core_remodeling",
                      "chr_te", "start_te", "end_te", "te", "len_te"] # rename columns

        df["enh_len"] = df["end_enh"] - df["start_enh"] # syn length

    else:
        if len(list(df)) ==16:
            df.columns = ["chr_syn", "start_syn", "end_syn", "enh_id",
                "chr_enh", "start_enh", "end_enh",
                          "seg_index", "core_remodeling", "core", "mrca",
                          "chr_te", "start_te", "end_te", "te", "len_te"] # rename columns


        else:
            df = df[[0,1,2,3,4,5,6,7,8,9,10,14,15,16,17,18]] # get rid of this mystery 4th column
            df.columns = ["chr_syn", "start_syn", "end_syn", "enh_id",
                "chr_enh", "start_enh", "end_enh",
                          "seg_index", "core_remodeling", "core", "mrca",
                          "chr_te", "start_te", "end_te", "te", "len_te"] # rename columns

        df["syn_len"] = df["end_syn"] - df["start_syn"] # syn length
        df["syn_id"] = df.chr_syn + ":" + df.end_syn.map(str) +"-"+ df.start_syn.map(str)# syn length

        df["enh_len"] = df["end_enh"] - df["start_enh"] # syn length

        df["code"] = "" # code the syntenic blocks
        df.loc[df["core_remodeling"] ==0, "code"] = "simple"
        df.loc[(df["core_remodeling"] ==1,) & (df.core ==1), "code"] = "complexcore"
        df.loc[(df["core_remodeling"] ==1,) & (df.core ==0), "code"] = "derived"

    df = df.loc[df.chr_enh != "chrX"]
    df = df.drop_duplicates() # drop duplicate
    df = df.loc[df["core_remodeling"] != 0.175] # idk what this is, but it is a BUG


    df = df.loc[df["enh_len"]>=6] # remove any syntenic blocks <6 bp long

    df = df.replace(".", 0) # replace all ".", with 0's
    df.len_te = df.len_te.fillna(0)
    df.len_te = df.len_te.astype(int) # change the datatype


    df["arch"] = "" # code the enhancers
    df.loc[df["core_remodeling"] ==0, "arch"] = "simple"
    df.loc[df["core_remodeling"] ==1, "arch"] = "complexenh"

    df["te_bin"] = 0 # te overlap binary

    df.loc[df["len_te"]>=6, "te_bin"] = 1 # re-write all the len_te that are less than 10bp long.
    df["dataset"] = data_id


    df.mrca=df.mrca.round(3)
    df = pd.merge(df, syn_gen_bkgd, how = "left", on = "mrca")
    return df


# In[6]:


def cleandf(outf):
    outpath = "/".join(outf.split("/")[:-1])

    outid = (outf.split("/")[-1]).split(".")[0]
    newf = "%scleaned_%s.bed"% (outpath, outid)
    cmd = '''awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' %s > %s''' % (outf, newf)

    subprocess.call(cmd, shell=True)
    return newf


# In[7]:


def fdr_formatting(results_dict):
    df = pd.concat(results_dict.values()) # concat results dict

    fdr_bool, fdr_pval =  multitest.fdrcorrection(df.pvalue, alpha=0.05, method='indep', is_sorted=False) # FDR 5%

    df["fdr_pval"] = fdr_pval # add FDR to df
    df["-log10p"] = np.log10(df["fdr_pval"])*-1 # -log10(p)

    df["log2_OR"]= np.log2(df["OR"]) # log2(OR)

    # count the number of TEs from family
    df["total_complex"] = df["complex_te"]+ df["complex_no_te"]
    df["total_simple"] = df["simple_te"]+ df["simple_no_te"]

    # plotting annotation
    df["fam_count_complex"] = df["fam"] + " (" + df["complex_te"].map(str) + "/"         + df["total_complex"].map(str) + ")"
    df["fam_count_simple"] = df["fam"] + " (" + df["simple_te"].map(str) + "/"         + df["total_simple"].map(str) + ")"

    return df


# In[8]:


in_df = pd.read_csv(outf, sep = '\t', header = None)
in_df.head()


# In[9]:


newf = cleandf(outf)
in_df = pd.read_csv(newf, sep = '\t', header = None)
in_df.head()


# In[10]:


len(in_df[3].unique())


# In[11]:


# Format enhdf
# add TE fam column
# make core_remodeling a float
# print # of unique enhancers overlapping TEs

enhdf = format_df(in_df, "FANTOM", TRIMMED)
enhdf["fam"]= ""
enhdf["fam"] = enhdf.loc[enhdf.te_bin!=0, "te"].apply(lambda x: "_".join(x.split("_")[1:]))
enhdf["fam"] = enhdf.loc[enhdf.te_bin!=0, "te"].apply(lambda x: (x.split("_")[-1]))
enhdf.core_remodeling=enhdf.core_remodeling.astype(float)
print(len(enhdf.enh_id.unique()))
enhdf.head()


# In[12]:


# groupby enhancer id, arch, te overlap
# count arch observations overlapping TEs
# FET simple v. complex enhancer TE overlap

enh = enhdf.groupby(["enh_id", "core_remodeling"])["te_bin"].max().reset_index()

obs = enh.groupby(["core_remodeling", "te_bin"])["enh_id"].count().reset_index()
print(obs)

simple_no_te = obs.enh_id.loc[(obs.core_remodeling == 0) &(obs.te_bin ==0)].item()
simple_te = obs.enh_id.loc[(obs.core_remodeling == 0) &(obs.te_bin ==1)].item()
complex_no_te = obs.enh_id.loc[(obs.core_remodeling == 1) &(obs.te_bin ==0)].item()
complex_te = obs.enh_id.loc[(obs.core_remodeling == 1) &(obs.te_bin ==1)].item()
fet_obs = [[ simple_te, simple_no_te], [complex_te, complex_no_te,]]

o, pvalue = stats.fisher_exact(fet_obs)
print(o, pvalue, fet_obs)


# In[13]:


RECONSTRUCT_SHUF_DICT = 0
shuf_dict={}
if TRIMMED ==1:

    shuf_te = glob.glob("%strimmed-310_shuf-*_x_te.bed"% outpath)

else:
    shuf_te = glob.glob("%sshuf-*.bed"% outpath)
if RECONSTRUCT_SHUF_DICT ==1:
    for shuf in shuf_te:

        shuf_id = (shuf.split("/")[-1]).split(".")[0]
        print(shuf_id)
        if "2501x" not in shuf_id:

            newf = cleandf(shuf)

            if os.path.getsize(shuf)>0:
                shuf_df = pd.read_csv(newf, sep = '\t', header = None, error_bad_lines = False)

                shufenh = format_df(shuf_df, shuf_id, 1)
                shufenh["fam"]= ""
                shufenh["fam"] = shufenh.te.loc[shufenh.te_bin!=0].apply(lambda x: "_".join(x.split("_")[1:]))
                shufenh["fam"] = shufenh.te.loc[shufenh.te_bin!=0].apply(lambda x: (x.split("_")[1:]))
                shuf_dict[shuf_id] = shufenh

        else:
            print("skipping", shuf_id)

        # enhancer level TE overlap - Weighted count, percent of length

    shufdf = pd.concat(shuf_dict.values())
    shufdf_outf = "%sall_trimmed-310_shuf.bed" % outpath
    shufdf.to_csv(shufdf_outf, sep ='\t', index = False)
else:
    shufdf_outf = "%sall_trimmed-310_shuf.bed" % outpath
    shufdf = pd.read_csv(shufdf_outf, sep ='\t')


# # get teenh

# In[14]:


def get_te_enh(df, dataset):
    teenh = df.loc[df.dataset == dataset].groupby(["enh_id", "core_remodeling", "dataset", "enh_len"])[["te_bin", "len_te"]].sum().reset_index()

    teenh["w_te_sum"] = teenh["te_bin"].divide(teenh["enh_len"]).round(4) # weighted # TE per length (i.e. how many TEs can you detect)

    teenh["per_te"] = teenh["len_te"].divide(teenh["enh_len"]).round(3) # % of length enh is TE?

    return teenh


# In[15]:


dsets = shufdf.dataset.unique()
dataset = dsets[0]
shufdf.loc[shufdf.dataset == dataset].groupby(["enh_id", "core_remodeling", "dataset", "enh_len"])[["te_bin", "len_te"]].sum().reset_index().head()


# In[16]:


shufdf.loc[shufdf.dataset == dataset].groupby(["enh_id", "core_remodeling"])[["te_bin", "len_te"]].sum().reset_index().head()


# In[17]:


teenh = get_te_enh(enhdf, "FANTOM")
teenh.head()

shufenh ={}
for i in shufdf.dataset.unique():
    print(i)
    test = get_te_enh(shufdf,i)
    shufenh[i] = test

plot = pd.concat(shufenh.values())
plot = pd.concat([plot, teenh])
plt.rcParams.update({'font.size': 15})
f, ax = plt.subplots(figsize = (16,8))
sns.pointplot(y = "w_te_sum",  x = "dataset", data = plot, hue = "core_remodeling",
           palette = es_pal, join = False)
ax.set_xlabel("")
ax.set_ylabel("# TEs per basepair")
ax.set_title("# TEs per FANTOM enh length")
ax.set_xticklabels(ax.get_xticklabels(), rotation = 270)
m, mp = stats.mannwhitneyu(teenh.w_te_sum.loc[teenh.core_remodeling ==0],
                           teenh.w_te_sum.loc[teenh.core_remodeling ==1])
ax.get_legend().remove()
print(m, mp)
plt.savefig('%sfantom_v_shuffle_te_arch_density.pdf' % RE, bbox_inches = "tight")

teenh.head()


# In[18]:


plot = pd.concat(shufenh.values())
plot.dataset.unique()


# In[19]:


bins = teenh.te_bin.loc[teenh.core_remodeling ==1].max()

sns.distplot(teenh.te_bin.loc[teenh.core_remodeling ==1], label = "complex", kde = False, bins =bins, norm_hist = True)
sns.distplot(teenh.te_bin.loc[teenh.core_remodeling ==0], label = "simple", kde = False, bins =bins,  norm_hist = True)
plt.legend()
print(teenh.groupby("core_remodeling")["te_bin"].mean())


# In[20]:


plot = pd.concat(shufenh.values())
plot = pd.concat([plot, teenh])
plt.rcParams.update({'font.size': 15})
f, ax = plt.subplots(figsize = (8,8))
sns.barplot(y = "w_te_sum",  x = "core_remodeling", data = plot.loc[plot.dataset.str.contains("FANTOM")],
           palette = es_pal,)
ax.set_xlabel("")
ax.set_ylabel("# TEs per basepair")
ax.set_title("# TEs per FANTOM enh length")
ax.set_xticklabels(ax.get_xticklabels(), rotation = 270)
m, mp = stats.mannwhitneyu(teenh.w_te_sum.loc[teenh.core_remodeling ==0],
                           teenh.w_te_sum.loc[teenh.core_remodeling ==1])
print(m, mp)
plt.savefig('%sfantom_te_arch_density_310fkb.pdf' % RE, bbox_inches = "tight")

teenh.groupby("core_remodeling")["w_te_sum"].mean()


# In[23]:


plt.rcParams.update({'font.size': 15})
f, ax = plt.subplots(figsize = (16,8))
sns.barplot(y = "per_te",  x = "dataset", data = plot, hue = "core_remodeling",
           palette = es_pal,)
ax.set_xlabel("")
ax.set_ylabel("TE as % enh")
ax.set_title("TEs as % FANTOM enh length")
m, mp = stats.mannwhitneyu(teenh.per_te.loc[teenh.core_remodeling ==0],
                           teenh.per_te.loc[teenh.core_remodeling ==1])
print(m, mp)
ax.set_xticklabels(ax.get_xticklabels(), rotation = 270)
plt.savefig('%sfantom_te_arch_percent_te_310.pdf' % RE, bbox_inches = "tight")

plot.head()


# In[24]:


plt.rcParams.update({'font.size': 15})
f, ax = plt.subplots(figsize = (8,8))
sns.barplot(y = "per_te",  x = "core_remodeling", data = plot.loc[plot.dataset.str.contains("FANTOM")],
           palette = es_pal,)
ax.set_xlabel("")
ax.set_ylabel("TE as % enh")
ax.set_title("% TE per FANTOM enh length")
m, mp = stats.mannwhitneyu(teenh.per_te.loc[teenh.core_remodeling ==0],
                           teenh.per_te.loc[teenh.core_remodeling ==1])
print(m, mp)
ax.set_xticklabels(ax.get_xticklabels(), rotation = 270)
plt.savefig('%sfantom_only_arch_percent_te_310.pdf' % RE, bbox_inches = "tight")

print(teenh.groupby("core_remodeling")["per_te"].mean())
plot.head()


# In[21]:


df = pd.concat([enhdf, shufdf])
df.head()


# In[22]:


pre_obs = df.groupby(["core_remodeling", "dataset", "enh_id"])["te_bin"].max().reset_index()
pre_obs["datatype"] = pre_obs.dataset.apply(lambda x: (x.split("-")[0]))

preobs = pre_obs.groupby(["core_remodeling", "datatype", "te_bin"])["enh_id"].count().reset_index()
preobs[["datatype", "core_remodeling", "te_bin", "enh_id"]].drop_duplicates()
preobs.head()


# # odds of complex, simple v. shuffles

# In[23]:


odds_dict ={}
# complex v background
obsf = preobs.loc[(preobs.datatype == "FANTOM") & (preobs.core_remodeling == 1) ]
obss = preobs.loc[(preobs.datatype.str.contains("trimmed"))& (preobs.core_remodeling == 1)]

fan_no_te = obsf.enh_id.loc[(obsf.te_bin ==0)].item()
fan_te = obsf.enh_id.loc[(obsf.te_bin ==1)].item()
shuf_no_te = obss.enh_id.loc[(obss.te_bin ==0)].item()
shuf_te = obss.enh_id.loc[(obss.te_bin ==1)].item()
fet_obs = [[ fan_te, fan_no_te], [shuf_te, shuf_no_te,]]

o, pvalue = stats.fisher_exact(fet_obs)
table = sm.stats.Table2x2(fet_obs) # get confidence interval
odds_ci = table.oddsratio_confint()

newdf = pd.DataFrame({"odds":[o],"p":[pvalue], "arch":["complex"],
                      "ci_lower": [odds_ci[0]],"ci_upper": [odds_ci[1]],})
odds_dict["complex"] = newdf

print("complex", o, pvalue)
print(fet_obs)

# complex v background
obsf = preobs.loc[(preobs.datatype == "FANTOM") & (preobs.core_remodeling == 0) ]
obss = preobs.loc[(preobs.datatype.str.contains("trimmed"))& (preobs.core_remodeling == 0) ]

fan_no_te = obsf.enh_id.loc[(obsf.te_bin ==0)].item()
fan_te = obsf.enh_id.loc[(obsf.te_bin ==1)].item()
shuf_no_te = obss.enh_id.loc[(obss.te_bin ==0)].item()
shuf_te = obss.enh_id.loc[(obss.te_bin ==1)].item()
fet_obs = [[ fan_te, fan_no_te], [shuf_te, shuf_no_te,]]

o, pvalue = stats.fisher_exact(fet_obs)
table = sm.stats.Table2x2(fet_obs) # get confidence interval
odds_ci = table.oddsratio_confint()

newdf = pd.DataFrame({"odds":[o],"p":[pvalue], "arch":["simple"],
                                            "ci_lower": [odds_ci[0]],"ci_upper": [odds_ci[1]],})
odds_dict["simple"] = newdf

print("simple", o, pvalue)
print(fet_obs)


# In[22]:


# core v. derived was calculated in fig5a-TE_x_FANTOM-310_syn_v_te_mrca.ipynb
o, pvalue, obs = 0.5642937173776406, 9.723858089297309e-89, [[2888, 8772], [4255, 7293]]
ci_lower, ci_upper = 0.5332697915946192, 0.5971225156401473
newdf = pd.DataFrame({"odds":[o],"p":[pvalue], "arch":["core"],
                      "ci_lower": [odds_ci[0]],"ci_upper": [odds_ci[1]],})
odds_dict["core"] = newdf


# # Figure 5a

# In[69]:


odds = pd.concat(odds_dict.values())
odds["log2"] = np.log2(odds.odds)

odds.head()


# In[70]:


odds.head()


# In[81]:


from matplotlib.ticker import MultipleLocator
import matplotlib.ticker as ticker

odds = pd.concat(odds_dict.values())
odds["log2"] = np.log2(odds.odds)

sns.set("poster")
fig, ax = plt.subplots(figsize = (8,8))

sns.barplot(x = "arch", y = "log2", data = odds,
            linewidth=2.5, facecolor=(1, 1, 1, 0),
             edgecolor=".2",  yerr=(odds["ci_upper"] - odds["ci_lower"]))
plt.axhline(0, color = "grey", linewidth = 2.5)

ax.set_ylabel("transposable element fold enrichment")
ax.set_xticklabels(["Complex\nEnhancer", "Simple\nEnhancer", "Complex\nEnhancer\nCore v. Derived"])

ax.set_xlabel("")

#ax.get_yaxis().ticker.LogLocator(base=2)

from matplotlib.ticker import MultipleLocator
import matplotlib.ticker as ticker
ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(2**x))
ax.yaxis.set_major_formatter(ticks)
ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.set_ylim(-2,0.5)
plt.savefig("%sfig5a.pdf" % RE, bbox_inches = "tight")


# In[62]:
