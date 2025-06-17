#!/usr/bin/env python
# coding: utf-8

# In[1]:


import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os, sys
import pybedtools as pb
from scipy import stats
import seaborn as sns

from statsmodels.stats import multitest
get_ipython().run_line_magic('matplotlib', 'inline')

colors = ["windows blue", "amber", "dusty purple", "faded green", "greyish"]
palette = sns.xkcd_palette(colors)
sns.palplot(palette)

es = ["faded green", "greyish"]
esPal = sns.xkcd_palette(es)

phast = ["light blue", "greyish"]
phastPal = sns.xkcd_palette(phast)

arch = ["dusty green", "amber"]
archPal = sns.xkcd_palette(arch)

plt.rcParams.update({'font.size': 18})


import datetime
LAST_RUN = datetime.datetime.now()
TODAY = (datetime.date.today())

freq = 0.1
RE = "/dors/capra_lab/projects/enhancer_ages/fantom/results/phastcons/"
if os.path.exists(RE) == False:
    os.system("mkdir %s" % RE)

print("last run", LAST_RUN)


# In[2]:


RUN_BED = 0


# In[3]:


path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/phastcons/"


# # load the tissue descriptions for FANTOM

# In[4]:


desc_file = "/dors/capra_lab/data/fantom/fantom5/facet_expressed_enhancers/sample_id_descriptions.txt"
desc_df= pd.read_csv(desc_file, sep = '\t', header = None)
desc_df.columns = ["sid", "desc"]
desc_dict = desc_df.to_dict("split")


# # syn background species

# In[5]:


syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t')
syn_gen_bkgd["mrca"]=syn_gen_bkgd["mrca"].round(3)
syn_gen_bkgd["mrca_2"]=syn_gen_bkgd["mrca_2"].round(3)
syn_gen_bkgd.head()


# # Phastcons bedfile

# In[6]:


# 46way vertebrate phastcons score.

phast_file = "/dors/capra_lab/data/evolutionary_conservation/phastcons/hg19/phastConsElements46wayVertebrate.bed"
p = pb.BedTool(phast_file).sort()


# # intersect Shuffles with Phastcons
# ### no exons, no x chromosomes, no enhancers longer than 10kb

# In[7]:


shuffle_path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/shuffle/breaks/"
shuffle_fs = glob.glob("%sshuf-all_fantom_enh*/*summary_matrix.bed" % shuffle_path)

shuf_dict = {}


for shuffle_file in shuffle_fs:

    shufval = int((shuffle_file.split("/")[-2]).split("enh-")[1])
    print(shufval)

    if shufval<100 and RUN_BED ==1: # evaluate 100 shuffles. Too much memory required for 1000

        shuffle_out_file = "%sshuffle-%s_enh_v_phastcons.bed" % (path, shufval)

        s = pb.BedTool(shuffle_file).sort() # turn shuf file into bed object

        s.intersect(p, f = freq, wao = True).saveas(shuffle_out_file) # intersect w/phastcons

        shuf_dict[shufval] = shuffle_out_file # add the outfile to the dictionary


# # intersect enhancers
# ### no exons, no x chromosome, no lengths longer than 10kb

# In[9]:


# TRIMMED FANTOM enhancers

file = "/dors/capra_lab/projects/enhancer_ages/fantom/data/all_fantom_enh/no_exon_all_fantom_summary_matrix.bed"

#e = pb.BedTool(file).sort()

fantom_phast_out_file = "%sfantomenh_v_phastcons.bed" % path # write the file

#e.intersect(p, f = freq, wao = True).saveas(fantom_phast_out_file)# intersection command


# # skip the Bedtools intersection

# In[10]:


shuf_dict.keys()


# In[15]:


fantom = pd.read_csv(fantom_phast_out_file, sep = '\t', header = None, usecols = [0,1,2,3,4,5,6,7,10,11,18,19])


# In[16]:


fantom.head()


# In[19]:


if RUN_BED ==0:
    fantom_phast_out_file = "%sfantomenh_v_phastcons.bed" % path
    shuffle_out_file = "%sshuffle_enh_v_phastcons.bed" % path

    # open dataframes

    cols = ["chr_enh", "start_enh", "end_enh",  "enh_id", "shuf_id",
             "seg_index", "core_remodeling", "arch", "mrca_2", "taxon2", "phastcon_score", "len_overlap"]

    fantom = pd.read_csv(fantom_phast_out_file, sep = '\t', header = None, usecols = [0,1,2,3,4,5,6,7,10,11,18,19])

    fantom.columns = cols
    fantom.shuf_id = "fantom"


    # gather all the shuffles
    shuf_dict = {}

    shuf_fs = glob.glob(f"{path}shuffle-*_enh_v_phastcons.bed")

    for f in shuf_fs:
        sid = (f.split("/")[-1]).split("_enh")[0]
        shuf_dict[sid] = f

    for key, val in shuf_dict.items():
        print(key, val)
        shuf_df = pd.read_csv(val, sep = '\t', header = None, usecols = [0,1,2,3,4,5,6,7,10,11,18,19])
        shuf_df.columns = cols
        shufdf[key] = shuf_df
    shuf_df = pd.concat(shufdf.values())


# In[20]:


#concatenate the shuffle dataframe
df = pd.concat([fantom, shuf_df])

### FILTER PHASTCONS MUST OVERLAP 10 BP OF ENHANCER ###
df["phastcon_overlap"]= np.where((df["len_overlap"]>=10),1,0)

df["data_type"] = "shuffle"
df.loc[df.shuf_id == "fantom", "data_type"] = "fantom"
df.head()


# In[21]:


xlabs = ["Homo", "Prim", "Euar", "Bore", "Euth",
        "Ther", "Mam", "Amni", "Tetr", "Vert"]
xlabs_nohomo = ["Prim", "Euar", "Bore", "Euth",
        "Ther", "Mam", "Amni", "Tetr", "Vert"]


# In[22]:


hue_order = ["fantom", "shuffle"]
sns.set("talk")
fig, ax = plt.subplots(figsize = (8,8))

sns.barplot(x = "taxon2", y = "phastcon_score",
              data = df.sort_values(by="mrca_2"),
              hue = "data_type", ax=ax, hue_order = hue_order,
              palette = esPal)

ax.set_xticklabels(xlabs, rotation = 90, horizontalalignment = "left")


# MWU age-distributions of fantom and shuffle

score_mwu, score_mwu_pval = stats.mannwhitneyu(df["phastcon_score"].loc[df["data_type"] == "fantom"],                            df["phastcon_score"].loc[df["data_type"] == "shuffle"])

ax.set(title="Phastcon-overlapping FANTOM Enhancers\nhave higher Phastcons_scores than 10x Genomic shuffle",
      xlabel = "Distribution Fantom/Shuffle PhastCons Scores\nMWU stat= %s, p-val = %s"% (score_mwu, score_mwu_pval))


plt.savefig("%sfantomenh_v_shuffle_phastcons_scores.pdf" % RE, bbox_inches = "tight")


# # Count of fantom and shuffle distributions overlapping phastCons

# In[23]:


len(df.loc[df.data_type == "fantom", "enh_id"].unique())


# In[24]:


code_grouped = df.groupby(["shuf_id", "mrca_2", "data_type", "taxon2","phastcon_overlap"])["enh_id"].count().reset_index()
print(code_grouped.head())

total_code = df.groupby(["shuf_id", "mrca_2"])["enh_id"].count().reset_index()
total_code.columns = [ "shuf_id", "mrca_2","total_code"]
total_code


# In[25]:


code_grouped = pd.merge(code_grouped, total_code, how = "left", on = (["shuf_id", "mrca_2"]))
code_grouped.head()


# In[26]:


code_grouped["freq"] = code_grouped["enh_id"].divide(code_grouped["total_code"])
code_grouped.loc[(code_grouped.shuf_id == "fantom") & (code_grouped.phastcon_overlap ==1)]


# In[32]:


shufg = code_grouped.loc[(code_grouped.shuf_id != "fantom") & (code_grouped.phastcon_overlap ==1)]
shufg.groupby(["mrca_2"])["enh_id"].sum()


# In[27]:


testStats = code_grouped.loc[code_grouped.phastcon_overlap ==1][["data_type", "mrca_2", "enh_id", "freq"]]

m, p = stats.mannwhitneyu(testStats.loc[(testStats.data_type == "fantom"), "freq"],
                  testStats.loc[(testStats.data_type.str.contains("shuf")), "freq"])

print(m, p)


# In[28]:


code_grouped.loc[(code_grouped.data_type == "fantom") & (code_grouped.phastcon_overlap ==1), "freq"].median()


# In[29]:


code_grouped.loc[(code_grouped.data_type.str.contains("shuf")) & (code_grouped.phastcon_overlap ==1), "freq"].median()


# In[104]:


enh = ["slate grey", "greyish"]
enhPal = sns.xkcd_palette(enh)

fig, ax = plt.subplots(figsize = (8,8))
sns.barplot(x = "taxon2",  y = "enh_id",
            data = code_grouped.loc[code_grouped.phastcon_overlap ==1].sort_values(by="mrca_2"),
            hue = "data_type",
           palette = enhPal)
ax.legend(bbox_to_anchor = (1,1))
ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
ax.set_ylabel("count phastCons overlap per age")
ax.set_xlabel("")

plt.savefig("%scounts_phastcons_fantom_v_shuffle.pdf" % RE, bbox_inches = "tight")


# In[105]:


enh = ["slate grey", "greyish"]
enhPal = sns.xkcd_palette(enh)

fig, ax = plt.subplots(figsize = (8,8))
sns.barplot(x = "taxon2",  y = "enh_id",
            data = code_grouped.loc[code_grouped.phastcon_overlap ==0].sort_values(by="mrca_2"),
            hue = "data_type",
           palette = enhPal)
ax.legend(bbox_to_anchor = (1,1))
ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
ax.set_ylabel("frequency within age\nNO phastCons overlap")
ax.set_xlabel("")

#plt.savefig("%sfig1e_counts_phastcons_fantom_v_shuffle.pdf" % RE, bbox_inches = "tight")


# In[106]:


enh = ["slate grey", "greyish"]
enhPal = sns.xkcd_palette(enh)

fig, ax = plt.subplots(figsize = (8,8))
sns.barplot(x = "taxon2",  y = "freq",
            data = code_grouped.loc[code_grouped.phastcon_overlap ==1].sort_values(by="mrca_2"),
            hue = "data_type",
           palette = enhPal)
ax.legend(bbox_to_anchor = (1,1))
ax.set_xticklabels(xlabs_nohomo, rotation = 90)
ax.set_ylabel("frequency within age\nphastCons overlap")
ax.set_xlabel("")

plt.savefig("%sfigureS1.2/FigS1.2e_phastcons_fantom_v_shuffle.pdf" % RE, bbox_inches = "tight")


# In[33]:


phastCorr = code_grouped.loc[code_grouped.phastcon_overlap ==1]
phastCorr


# In[ ]:





# In[108]:


xs = phastCorr.loc[phastCorr.data_type.str.contains("fantom"), "mrca_2"]
ys = phastCorr.loc[phastCorr.data_type.str.contains("fantom"), "freq"]

slope_s, intercept_s, r_value_s, p_value_s, std_err_s = stats.linregress(xs,ys)
print(slope_s, intercept_s, r_value_s, p_value_s, std_err_s)

xsShuf = phastCorr.loc[phastCorr.data_type.str.contains("shuffle"), "mrca_2"]
ysShuf = phastCorr.loc[phastCorr.data_type.str.contains("shuffle"), "freq"]

slope_s, intercept_s, r_value_s, p_value_s, std_err_s = stats.linregress(xsShuf,ysShuf)
print(slope_s, intercept_s, r_value_s, p_value_s, std_err_s)


# # How many fantom and shuffled enhancers (Simple | Complexenh) overlap PhastCons Elements?
#
# ## Count
# ## Frequency

# In[34]:


fantom = df.loc[(df.data_type == "fantom") & (df.core_remodeling != 0.175)]

fantom.head()


# In[35]:


len(fantom.enh_id.unique())


# In[36]:


fantom.groupby(["core_remodeling", "mrca_2", ])["enh_id"].count()


# In[37]:


colors = [ "amber", "dusty green",]
pal = sns.xkcd_palette(colors)


# In[38]:


fsum = fantom.groupby("core_remodeling")["phastcon_overlap"].sum().reset_index()
fsum.columns = ["core_remodeling", "phastSum"]

fcount = fantom.groupby("core_remodeling")["phastcon_overlap"].count().reset_index()
fcount.columns = ["core_remodeling", "phastCounts"]

fresults = pd.merge(fsum, fcount)
fresults["freq"] = fresults.phastSum.divide(fresults.phastCounts)
fresults


# In[39]:


labs = ["Simple", "Complex"]
order = [0,1]
fig, (ax) = plt.subplots( figsize = (8,8))
sns.set("poster")
sns.set_style("white")
sns.barplot(x = "core_remodeling", y = "freq",
            data = fresults, palette = pal,
           order = order)

ax.set_xticklabels(labs)
ax.set(xlabel = "",
       ylabel = "% architecture\nin phastcon")
plt.savefig("%sfigures3.1/FigS3.1fantomenh_arch_phast_percent.pdf"%RE, bbox_inches = "tight")


fig, (ax) = plt.subplots( figsize = (8,8))
sns.set("poster")
sns.set_style("white")
sns.barplot(x = "core_remodeling", y = "phastSum",
            data = fresults, palette = pal,
           order = order)

ax.set_xticklabels(labs)
ax.set(xlabel = "",
       ylabel = "% architecture\nin phastcon")
plt.savefig("%sfigures3.1/FigS3.1fantomenh_arch_phast_count.pdf"%RE, bbox_inches = "tight")


# In[40]:


fsum = fantom.groupby(["core_remodeling", "mrca_2", "taxon2"])["phastcon_overlap"].sum().reset_index()
fsum.columns = ["core_remodeling", "mrca_2", "taxon2", "phastSum"]

fcount = fantom.groupby(["core_remodeling", "mrca_2", "taxon2"])["phastcon_overlap"].count().reset_index()
fcount.columns = ["core_remodeling", "mrca_2", "taxon2", "phastCounts"]

fresults = pd.merge(fsum, fcount)
fresults["freq"] = fresults.phastSum.divide(fresults.phastCounts)


# In[44]:


xlabs = ['Prim', 'Euar', 'Bore', 'Euth', 'Ther', 'Mam', 'Amni', 'Tetr', 'Vert']


# In[45]:


fig, (ax1) = plt.subplots( figsize = (8,8))
sns.set("poster")
sns.set_style("white")
sns.barplot(x = "taxon2", y = "freq",
            data = fresults.sort_values(by="mrca_2"),
            hue = 'core_remodeling', palette = pal, ax = ax1)

ax1.set_xticklabels(xlabs , rotation = 90)

ax1.set(ylabel = "% Phastcons overlap per age", xlabel = "")

ax1.legend().remove()

plt.savefig("%sfigS10fantomenh_arch_phast_mrca_x_percent.pdf"%RE, bbox_inches = "tight")


fig, (ax2) = plt.subplots( figsize = (8,8))
sns.barplot(x = "taxon2", y = "phastSum",
            data = fresults.sort_values(by="mrca_2"), hue = 'core_remodeling', palette = pal, ax=ax2)
ax2.set(ylabel= "Count Phastcons overlap", xlabel = "")
ax2.legend().remove()
ax2.set_xticklabels(xlabs, rotation = 90)

plt.savefig("%sfigS10_fantomenh_arch_phast_count.pdf"%RE, bbox_inches = "tight")


# In[ ]:
