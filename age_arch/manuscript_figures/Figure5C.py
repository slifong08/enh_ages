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
import subprocess

today = datetime.date.today()
RE = "/dors/capra_lab/projects/enhancer_ages/fantom/results/te/"
if os.path.exists(RE) == False:
    cmd = "mkdir %s" %RE
    os.system(cmd)

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


simp = ["lemon", "amber", "mustard", "greyish"]
simp_pal = sns.xkcd_palette(simp)
sns.palplot(simp_pal)

comp = ["pale green", "faded green", "green grey", "greyish"]
comp_pal = sns.xkcd_palette(comp)
sns.palplot(comp_pal)
arch_te = [ "greyish", "sky blue"]
arch_te_pal = sns.xkcd_palette(arch_te)
sns.palplot(arch_te_pal)

arch_te = [ "greyish", "sky blue"]
arch_te_pal = sns.xkcd_palette(arch_te)
sns.palplot(arch_te_pal)

simp2 = ["amber", "greyish"]
simp2_pal = sns.xkcd_palette(simp2)
sns.palplot(simp2_pal)

comp2 = ["faded green",  "greyish"]
comp2_pal = sns.xkcd_palette(comp2)
sns.palplot(comp2_pal)


# In[3]:


TRIMMED = 310 # run syntenic block file or enhancer file?
BED_INTERSECTION = 0

te_path = "/dors/capra_lab/users/abraha1/projects/transposable_elements/data/"
tefam = "%shg19_TE_counts_wlin.txt" % te_path
dfam_file="%sdfam/hg38_dfam.nrph.hits_genomicCoordinates_tidy.bed" % te_path
repeatmasker = "%sfiltered_formated_hg19fromhg38-TE_coords.tsv" % te_path

syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2", "mya", "mya2"]] # whittle down the df

# handle the TEs that have an old age
old = pd.DataFrame({"mrca": [0.39], "taxon": ["older_than_mammalia"], "mrca_2": [0.39],
                        "taxon2": ["older_than_mammalia"], "mya": [">177"], "mya2": [">177"]})

syn_gen_bkgd = pd.concat([syn_gen_bkgd, old])


# ### FORMAT TE FAMILY FILE
# ### to match syntenic background

# In[4]:


famdf = pd.read_csv(tefam, sep = '\t')
famdf.columns = ["te_fam", "fam", "taxon", "count"]
famdf.taxon.loc[famdf.taxon == "Primates"] = "Primate" # make consisten w/ syn_gen_bkgd
famdf.taxon.loc[famdf.taxon == "theria"] = "Theria" # make consisten w/ syn_gen_bkgd
famdf.taxon.loc[famdf.taxon == "Euarchonta"] = "Euarchontoglires"
famdf.taxon.loc[famdf.taxon == "Hominoidea"] = "Hominidae"
famdf.taxon.loc[famdf.taxon == "Rodentia"] = "Euarchontoglires"
famdf.taxon.loc[famdf.taxon == "Muridae"] = "Euarchontoglires"
famdf.taxon.loc[famdf.taxon == "Homo_sapiens"] = "Homo sapiens"# make consisten w/ syn_gen_bkgd
famdf.taxon.loc[famdf.taxon == "old"] = "older_than_mammalia"

famdf = pd.merge(famdf, syn_gen_bkgd[["mrca", "taxon", "mrca_2",  "taxon2"]], how = "left", on = "taxon") # get ages


famdf.columns = ["te", "fam", "taxon_te", "count_te",  "mrca_te",  "mrca_2_te",  "taxon2_te"]
famdf.head()


# In[5]:


famdf.mrca_te.unique()


# In[6]:


# INTERSECT FANTOM
fantom_path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/"
outpath= "/dors/capra_lab/projects/enhancer_ages/fantom/data/te/"

fantom_path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/te/"
outpath = fantom_path

if TRIMMED ==310:
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


# In[7]:


# shuffle x TE intersection

BED_INTERSECTION = 0

if TRIMMED == 310:
    shuf_path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/shuffle/breaks/"

    shuffles = glob.glob("%strimmed-310_shuf-*.bed"% shuf_path)

elif TRIMMED == 1000:

    shuf_path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/shuffle/breaks/"

    shuffles = glob.glob("%strimmed-1000_shuf-*.bed"% shuf_path)

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
            outf = "%s%s_x_te.bed" % (shuf_path, shuf_id) # DORS
            outf = "%s%s_x_te.bed" % (homepath, shuf_id) # Home

            s.intersect(r, wao = True).saveas(outf)


# In[8]:


def fdr_formatting(results_dict):
    df = pd.concat(results_dict.values()) # concat results dict

    fdr_bool, fdr_pval =  multitest.fdrcorrection(df.pvalue, alpha=0.05, method='indep', is_sorted=False) # FDR 5%

    df["fdr_pval"] = fdr_pval # add FDR to df
    df["-log10p"] = np.log10(df["fdr_pval"])*-1 # -log10(p)

    df["log2_odds"]= np.log2(df["odds"]) # log2(OR)

    return df


# In[9]:


def format_df(df, data_id, trimmed):
    if trimmed >0:
        df.columns = ["chr_enh", "start_enh", "end_enh", "old_enh_len", "old_enh_id", "mrca", "core_remodeling",
                     "chr_te", "start_te", "end_te", "te_fam", "len_te"]


    else:
        df = df[[0,1,2,3,4,5,6,7,8,9,10,14,15,16,17,18]] # get rid of this mystery 4th column
        df.columns = ["chr_syn", "start_syn", "end_syn", "enh_id",
            "chr_enh", "start_enh", "end_enh",
                      "seg_index", "core_remodeling", "core", "mrca",
                      "chr_te", "start_te", "end_te", "te_fam", "len_te"] # rename columns

        df["syn_len"] = df["end_syn"] - df["start_syn"] # syn length
        df["syn_id"] = df.chr_syn + ":" + df.end_syn.map(str) +"-"+ df.start_syn.map(str)# syn length
        df = df.loc[df["syn_len"]>=6] # remove any syntenic blocks <6 bp long


        df["code"] = "" # code the syntenic blocks
        df["code"].loc[df["core_remodeling"] ==0] = "simple"
        df["code"].loc[(df["core_remodeling"] ==1) & (df.core ==1)] = "complexcore"
        df["code"].loc[(df["core_remodeling"] ==1) & (df.core ==0)] = "derived"

    df["enh_id"] = df.chr_enh + ":" + df.end_enh.map(str) +"-"+ df.start_enh.map(str)# syn length
    df["enh_len"] = df["end_enh"] - df["start_enh"] # syn length
    df = df.loc[df.old_enh_len>=6]
    df = df.loc[df.chr_enh != "chrX"]
    df = df.drop_duplicates() # drop duplicate

    df = df.loc[df["core_remodeling"] != 0.175] # idk what this is, but it is a BUG


    df = df.replace(".", 0) # replace all ".", with 0's
    df.core_remodeling=df.core_remodeling.astype(float)
    df.len_te = df.len_te.fillna(0)
    df.len_te = df.len_te.astype(int) # change the datatype


    df["arch"] = "" # code the enhancers
    df["arch"].loc[df["core_remodeling"] ==0] = "simple"
    df["arch"].loc[df["core_remodeling"] ==1] = "complexenh"

    df["te_bin"] = 0 # te overlap binary

    df.te_bin.loc[df["len_te"]>= 20] = 1 # re-write all the len_te that are less than 10bp long.
    df["dataset"] = data_id


    df.mrca=df.mrca.round(3)

    df["fam"]= ""
    df["fam"] = df.te_fam.loc[df.te_fam != 0].apply(lambda x: (x.split("_")[-1]))
    df["te"] = df.te_fam.loc[df.te_fam != 0].apply(lambda x: "".join(x.split("_")[0]))
    df = pd.merge(df, syn_gen_bkgd, how = "left", on = "mrca")

    return df


# In[10]:


def cleandf(outf):
    outpath = "/".join(outf.split("/")[:-1])

    outid = (outf.split("/")[-1]).split(".")[0]
    newf = "%scleaned_%s.bed"% (outpath, outid)
    cmd = '''awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' %s > %s''' % (outf, newf)

    subprocess.call(cmd, shell=True)
    return newf


# In[11]:


newf = cleandf(outf)
in_df = pd.read_csv(newf, sep = '\t', header = None)
print(len(list(in_df)))

# Format enhdf

enhdf = format_df(in_df, "FANTOM", TRIMMED)

# get oldest enhancer age
enh_mrca = enhdf.groupby(["enh_id"])["mrca", "mrca_2"].max().reset_index()

enh_mrca = pd.merge(enh_mrca, syn_gen_bkgd, how = "left")

enh_mrca.columns = ["enh_id", "core_mrca", "core_mrca_2", "mya", "mya2", "taxon", "taxon2",]

# merge enhancer+TE intersection with TE family information

enh_te = pd.merge(enhdf[["enh_id", "old_enh_id", "te", "fam", "te_fam","core_remodeling",
                             "te_bin", "len_te", "mrca", "mrca_2"]],
                    famdf, how = "left")\
                    .drop_duplicates()

# merge enh+TE information with oldest enhancer ages information
enh_te = pd.merge (enh_te, enh_mrca, how = "left", on = "enh_id")
enh_te.head()


# In[12]:


enh_te["class"] ="No TE"
enh_te["class"].loc[(enh_te.core_mrca == enh_te.mrca) &(enh_te.te_bin ==1)] = "TE overlap"
enh_te["class"].loc[(enh_te.core_mrca<enh_te.mrca)  &(enh_te.te_bin ==1)] = "TE older"
enh_te["class"].loc[(enh_te.core_mrca>enh_te.mrca)  &(enh_te.te_bin ==1)] = "TE younger"

enh_te["class2"] ="No TE"
enh_te["class2"].loc[(enh_te.core_mrca_2 == enh_te.mrca_2) &(enh_te.te_bin ==1)] = "TE overlap"
enh_te["class2"].loc[(enh_te.core_mrca_2<enh_te.mrca_2)  &(enh_te.te_bin ==1)] = "TE older"
enh_te["class2"].loc[(enh_te.core_mrca_2>enh_te.mrca_2)  &(enh_te.te_bin ==1)] = "TE younger"

enh_te.head()


# In[13]:


oldest = enh_te.groupby(["core_remodeling", "enh_id"])[["core_mrca", "mrca", "mrca_te", "te_bin"]].max().reset_index()
oldest["class"] ="No TE"
oldest["class"].loc[(oldest.core_mrca == oldest.mrca) &(oldest.te_bin ==1)] = "TE overlap"
oldest["class"].loc[(oldest.core_mrca<oldest.mrca) &(oldest.te_bin ==1)] = "TE older"
oldest["class"].loc[(oldest.core_mrca>oldest.mrca) &(oldest.te_bin ==1)] = "TE younger"

plotenh = oldest.groupby(["core_remodeling", "mrca", "class"])["enh_id"].count().reset_index()
enh_totals = plotenh.groupby("core_remodeling")["enh_id"].sum().reset_index()
enh_totals.columns = ["core_remodeling", "core_remodeling_totals"]
plotenh = pd.merge(plotenh, enh_totals, how = "left")
plotenh["frac_arch"] = plotenh.enh_id.divide(plotenh.core_remodeling_totals)
plotenh = pd.merge(plotenh, syn_gen_bkgd, how = "left", on = "mrca")
plotenh.head()


# In[14]:


plotenh["frac_arch"].sum()


# In[18]:


sns.set_style("whitegrid")
sns.set(font_scale=2)
oldest = enh_te.groupby(["core_remodeling", "enh_id"])[["core_mrca_2", "mrca_2", "mrca_2_te", "te_bin"]].max().reset_index()

oldest["class"] ="No TE"
oldest["class"].loc[(oldest.core_mrca_2 == oldest.mrca_2) &(oldest.te_bin ==1)] = "TE overlap"
oldest["class"].loc[(oldest.core_mrca_2<oldest.mrca_2) &(oldest.te_bin ==1)] = "TE older"
oldest["class"].loc[(oldest.core_mrca_2>oldest.mrca_2) &(oldest.te_bin ==1)] = "TE younger"

plotenh = oldest.groupby(["core_remodeling", "mrca_2", "class"])["enh_id"].count().reset_index()

enh_totals = plotenh.groupby("core_remodeling")["enh_id"].sum().reset_index()

enh_totals.columns = ["core_remodeling", "core_remodeling_totals"]

plotenh = pd.merge(plotenh, enh_totals, how = "left")

plotenh["frac_arch"] = plotenh.enh_id.divide(plotenh.core_remodeling_totals)

plotenh = pd.merge(plotenh, syn_gen_bkgd[["mrca_2", "taxon2"]],
                       how = "left", on = "mrca_2").drop_duplicates()
plotenh.loc[plotenh["class"] == "TE overlap", "class"] = "TE overlap"
print(plotenh["frac_arch"].sum())

order = ["TE overlap","No TE",]


fig, (ax1, ax2) = plt.subplots(ncols = 2, figsize = (16,6))
plt.rcParams.update({'font.size': 15})


sns.barplot(x ="taxon2" , y = "frac_arch",
            data = plotenh.loc[plotenh["core_remodeling"]==1].sort_values(by = "mrca_2"),
            hue = "class", ax = ax1, hue_order = order,
            linewidth=2.5,  edgecolor=".2",
            palette = comp_pal)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation = 270)
ax1.set_ylabel("Complex Enhancer Frequency")
ax1.set_xlabel("oldest enhancer sequence age")
ax1.legend(loc = "upper right", frameon = False)

sns.barplot(x = "taxon2", y = "frac_arch",
            data = plotenh.loc[plotenh["core_remodeling"]==0].sort_values(by = "mrca_2"), ax =ax2,
            linewidth=2.5,  edgecolor=".2",
            hue = "class", hue_order = order, palette = simp_pal)

ax2.set_xticklabels(ax2.get_xticklabels(), rotation = 270)
ax2.set_ylabel("Simple Enhancer Frequency")
ax2.set_xlabel("oldest enhancer sequence age")
ax2.legend(loc = "upper right", frameon = False)
ax1.set_ylim(0,0.4)
ax2.set_ylim(0,0.4)
sns.set_style("white")
sns.set(font_scale=2)

plt.savefig("%sfigs5c-fantom_trimmed_310_mrca2_v_te_mrca_freq.pdf" % RE, bbox_inches = "tight")


# In[19]:


rel = plotenh.loc[plotenh["class"] == "TE overlap"]
rel.head()


# In[21]:


plotenh.head()


# In[25]:


age_totals = plotenh.groupby(["core_remodeling", "mrca_2", "taxon2"])["enh_id"].sum().reset_index()
age_totals.columns = ["core_remodeling", "mrca_2", "taxon2", "age_arch_totals"]


# In[28]:


age_totals = pd.merge(rel, age_totals)
age_totals["frac_age"] = age_totals.enh_id.divide(age_totals.age_arch_totals)
age_totals


# In[42]:
def get_counts(res, strat):

    if strat == 1:
        counts = res[["mrca_2", "core_remodeling", "enh_id"]].drop_duplicates()

        # add empty dataframe for complex human enhancers (do not exist by our def)
        empty = pd.DataFrame({"mrca_2": [0.000],
        "core_remodeling": [1],
        "enh_id": [0]
        })

        counts = pd.concat([counts, empty]) # concat the dataframe

        counts =counts.sort_values(by = ["core_remodeling", 'mrca_2']).reset_index()

    else:
        counts = res.groupby("arch")["enh_id"].count().reset_index()
        counts =counts.sort_values(by = "arch", ascending = False).reset_index()

    counts["enh_id"] = counts["enh_id"].astype(int)
    counts = counts.drop(["index"], axis = 1)

    return counts


fig, (ax) = plt.subplots(figsize = (6,6))
plt.rcParams.update({'font.size': 15})
sns.set_style("white")
sns.set(font_scale=2)

x = "mrca_2"
y = "frac_age"
data = age_totals
hue = "core_remodeling"

xlabs = ["Homo", "Prim", "Euar", "Bore", "Euth", "Ther", "Mam", "Amni", "Tetr", "Vert"]

splot = sns.barplot(x =x, y=y,
            data =data,
            hue = hue,
            ax = ax,
            linewidth=2.5,  edgecolor=".2",
            palette = es_pal)
STRAT = 1
counts = get_counts(age_totals, STRAT)
for n, p in enumerate(splot.patches):
    value = counts.iloc[n]["enh_id"].astype(int)
    splot.annotate(value,
                   (p.get_x() + p.get_width() / 2.,0.02),
                   ha = 'center', va = 'baseline',
                   size=15,
                   rotation = 90,
                   color = "white",
                   xytext = (0, 1),
                   textcoords = 'offset points'
                   )
ax.set_xticklabels(xlabs, rotation = 90)
ax.set(ylabel = "% Arch overlapping TE\nper age", xlabel = "")
ax.legend(bbox_to_anchor =(1,1), frameon = False).remove()


plt.savefig("%sFig5c-fantom_trimmed_310_mrca2_v_te_mrca_freq.pdf" % RE, bbox_inches = "tight")


# In[32]:


es_colors = ["amber", "faded green"]
palette = sns.xkcd_palette(es_colors)
order = [0.0 ,1]

fig, ax1 = plt.subplots(figsize = (6,6))
plt.rcParams.update({'font.size': 15})

sns.barplot(x ="taxon2" , y = "frac_arch",
            data = rel.sort_values(by = "mrca_2"),
            hue = "core_remodeling", ax = ax1, hue_order = order,
            linewidth=2.5,  edgecolor=".2",
            palette = palette)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation = 270)
ax1.set_ylabel("TE Enhancer Frequency")
ax1.set_xlabel("oldest enhancer sequence age")
ax1.legend(loc = "upper right", frameon = False)

ax1.set_ylim(0,0.4)

sns.set_style("white")
sns.set(font_scale=2)

plt.savefig("%sfigs5c-fantom_trimmed_310_mrca2_v_te_mrca_freq.pdf" % RE, bbox_inches = "tight")


# # MWU on COMPLEX enhancer MRCAs +/- TE

# In[23]:


statsComplex = oldest.loc[oldest.core_remodeling ==1]

m, p = stats.mannwhitneyu(statsComplex.loc[statsComplex.te_bin ==1, "core_mrca_2"],
                  statsComplex.loc[statsComplex.te_bin ==0, "core_mrca_2"],)
print(m,p)


# In[24]:


statsComplex.loc[statsComplex.te_bin ==1, "core_mrca_2"].median()


# In[25]:


statsComplex.loc[statsComplex.te_bin ==0, "core_mrca_2"].median()


# # MWU on SIMPLE enhancer MRCAs +/- TE

# In[26]:


# MWU on enhancer MRCAs +/- TE
statsSimple = oldest.loc[oldest.core_remodeling ==1]
m, p = stats.mannwhitneyu(statsSimple.loc[statsSimple.te_bin ==1, "core_mrca_2"],
                  statsSimple.loc[statsSimple.te_bin ==0, "core_mrca_2"],)
print(m,p)


# In[27]:


statsSimple.loc[statsSimple.te_bin ==1, "core_mrca_2"].median()


# In[28]:


statsSimple.loc[statsSimple.te_bin ==0, "core_mrca_2"].median()


# In[29]:


plotenh.head()
