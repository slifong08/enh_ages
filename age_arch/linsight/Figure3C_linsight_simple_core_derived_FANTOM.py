#!/usr/bin/env python
# coding: utf-8

# In[1]:


import glob
import pandas
import matplotlib.pyplot as plt
import numpy as np
import os, sys
from scipy import stats
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')

# sns colors
arch_colors = [ "amber", "dusty purple", "windows blue","greyish"]
arch_palette = sns.xkcd_palette(arch_colors)
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
RE = "/dors/capra_lab/projects/enhancer_ages/linsight/results/"


# In[2]:


linsight_path = "/dors/capra_lab/projects/enhancer_ages/linsight/data/"
fantom_fs = glob.glob("%sFANTOM_*.bed" % linsight_path)

# enh only
fantom_fs = "/dors/capra_lab/projects/enhancer_ages/linsight/data/all_unique_fantom_erna_112_tissue_linsight.bed"


# In[3]:


syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pandas.read_csv(syn_gen_bkgd_file, sep = '\t')
syn_gen_bkgd["mrca"]=syn_gen_bkgd["mrca"].round(3)
syn_gen_bkgd["mrca_2"]=syn_gen_bkgd["mrca_2"].round(3)
syn_gen_bkgd.head()
syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2"]]


# In[4]:



arch_id = (fantom_fs.split("/")[-1]).split("_")[1]
print(arch_id)
df = pandas.read_csv(fantom_fs, sep = '\t', header = None, low_memory=False)

df.columns = ['chr_syn','start_syn','end_syn', 'enh_id',
              'chr_enh', 'start_enh','end_enh',
              'seg_index', 'core_remodeling', 'core',
              'mrca', 'code', 'syn_id',
              "chr_lin",  "start_lin", "end_lin","linsight_score", "overlap"]


# In[5]:


df.head()


# In[6]:


df['enh_len'] = df.end_enh - df.start_enh
df['syn_len'] = df.end_syn - df.start_syn
df.head()


# In[7]:


# arch code
df["code"] = ""
df.code.loc[(df.core_remodeling == 0)& (df.core == 1)] = "simple"
df.code.loc[(df.core_remodeling == 1)& (df.core == 1)] = "complex_core"
df.code.loc[(df.core_remodeling == 1)& (df.core == 0)] = "derived"
df["arch"] = ""
df.arch.loc[(df.core_remodeling == 0)] = "simple"
df.arch.loc[(df.core_remodeling == 1)] = "complexenh"

df = df.loc[df.linsight_score != "."] # exclude the loci that do not overlap a linsight score

df.linsight_score = df.linsight_score.astype(float) # turn linsight scores into floats

df["linsight_id"] = df.chr_lin + ":" + df.start_lin.map(str) +"-"+ df.end_lin.map(str)

# make a dataframe only based on linsight scores and enhancer architectures

base_df = df[["chr_lin", "start_lin", "end_lin", "linsight_score", "code", "arch"]].drop_duplicates()

base_df["lin_len"] = base_df.end_lin - base_df.start_lin

base_df.head()

# expand linsight score estimates across bases

simple_df = base_df.loc[base_df.code.str.contains("simple")]
simple = np.repeat(simple_df.linsight_score, simple_df.lin_len) # expand linsight value for each simple basepair

derived_df = base_df.loc[base_df.code.str.contains("derived")]
derived = np.repeat(derived_df.linsight_score, derived_df.lin_len)

core_df = base_df.loc[base_df.code.str.contains("core")]
core = np.repeat(core_df.linsight_score, core_df.lin_len)

complexenh = pandas.concat([derived, core])
len(complexenh)


# In[8]:


len(df.enh_id.unique())


# In[9]:


df.shape


# In[10]:


df.groupby(["enh_id"])


# In[11]:


fig, ax = plt.subplots(figsize = (8, 8))

sns.distplot(derived, hist_kws=dict(cumulative=True),
                 kde_kws=dict(cumulative=True), label = "derived", kde =False, norm_hist = True, color = "blue")
sns.distplot(core, hist_kws=dict(cumulative=True),
                 kde_kws=dict(cumulative=True), label = "complex core", kde =False, norm_hist = True, color = "purple")
sns.distplot(simple, hist_kws=dict(cumulative=True),
                 kde_kws=dict(cumulative=True), label = "simple", kde =False, norm_hist = True, color = "gold")

k, kp = stats.kruskal(simple, derived, core)

ax.set_title("FANTOM LINSIGHT scores\nCumulative Distribution")
ax.set_xlabel("LINSIGHT score\nkruskal = %s, p = %s" % (k, kp))
ax.set_ylabel("% of enhancer bases")
ax.legend(bbox_to_anchor=(1.5, 1.0))

plt.savefig("%sfantom_linsight_architecture.pdf" %(RE), bbox_inches = "tight")


# In[12]:


simple


# In[13]:


simplef = simple.to_frame()
simplef["code"] = "simple"
simplef["arch"] = "simple"

derivedf = derived.to_frame()
derivedf["code"] = "derived"
derivedf["arch"] = "complex"

coref = core.to_frame()
coref["code"] = "complex_core"
coref["arch"] = "complex"

concat = pandas.concat([simplef, derivedf, coref])


# In[14]:


concat.groupby("code")["linsight_score"].median()


# In[15]:


concat.groupby("code")["linsight_score"].mean()


# In[16]:


order = ["simple", "complex_core", "derived"]
fig, ax = plt.subplots(figsize = (8, 8))
sns.boxplot(x = "code", y = "linsight_score", data = concat,
            showfliers = False, order = order, notch = True,
           palette = arch_palette)
ax.set_xlabel("\nkruskal = %s, p = %s" % (k, kp))
#ax.set_title("LINSIGHT score by FANTOM architecture age")
ax.set_xlabel("")
ax.set_ylabel("LINSIGHT score")
sns.set("poster")
plt.savefig("%sfantom_linsight_architecture_boxplot.pdf" %(RE), bbox_inches = "tight")


# In[17]:


colores = ["amber", "faded green"]
pal = sns.xkcd_palette(colores)
sample = concat.sample(frac =0.03)
order = ["simple", "complex",]
plt.subplot(1,2,2)
sns.boxplot(x = "arch", y = "linsight_score", data = concat,
            showfliers = False, order = order, notch = True,
           palette = pal)

ax.set_xlabel("")
ax.set_xticklabels("")
ax.set_ylabel("LINSIGHT score")
ax.set_ylim(0,0.65)
sns.set("poster")
plt.subplot(2,2,2)
sample_points = base_df.sample(frac =0.03)
sns.pointplot(x = "arch", y = "linsight_score", data = sample_points,
              hue = "arch",
              #palette = pal,
              hue_order = order, join = False, dodge = 0.25)
sns.set("poster")
ax.set_xticklabels(ax.get_xticklabels(), rotation = 270, horizontalalignment = "left")
ax.set_xlabel("architecture age")

ax.legend().remove()
ax.set_xlabel("")
ax.set_xticklabels("")
ax.set_ylabel("LINSIGHT score")


# In[ ]:


sample.groupby("arch").count()


# In[ ]:


concat.groupby("arch")["linsight_score"].median()


# In[ ]:


concat.groupby("arch")["linsight_score"].mean()


# In[ ]:



# make a dataframe only based on linsight scores and enhancer architectures

base_df = df[["linsight_id", "enh_id" ,"chr_lin", "start_lin", "end_lin", "linsight_score", "mrca", "code", "arch"]].drop_duplicates()
base_df.mrca = base_df.mrca.round(3)
base_df["lin_len"] = base_df.end_lin - base_df.start_lin

base_df = base_df.loc[base_df.code != ""] # exclude non-coded regions

base_df = pandas.merge(base_df, syn_gen_bkgd) # add in taxon2

base_df.sort_values(by="mrca_2").head()
base_df.mrca_2.unique()


# In[ ]:


core_age = base_df.groupby("enh_id")["mrca_2"].max().reset_index()
core_age.columns = ["enh_id", "core_mrca"]
base_df = pandas.merge(base_df, core_age, how = 'left', on = 'enh_id')
base_df.head()


# In[ ]:


order = ["simple", "complex_core", "derived"]
fig, ax = plt.subplots(figsize = (8, 8))
sns.barplot(x = "taxon2", y = "linsight_score", data = base_df.sort_values(by="mrca_2")
              , hue = "code", palette = arch_palette, hue_order = order)
ax.set_xticklabels(ax.get_xticklabels(), rotation = 270, horizontalalignment = "left")
ax.set_xlabel("Sequence age")
ax.set_title("LINSIGHT score by FANTOM architecture age")
ax.legend(bbox_to_anchor=(1.5, 1.0))
plt.savefig("%sfantom_linsight_architecture_mrca.pdf" %(RE), bbox_inches = "tight")


# In[ ]:


order = ["simple", "complex_core", "derived"]
fig, ax = plt.subplots(figsize = (8, 8))
sns.barplot(x = "core_mrca", y = "linsight_score", data = base_df.sort_values(by="mrca_2")
              , hue = "code", palette = arch_palette, hue_order = order)
ax.set_xticklabels(ax.get_xticklabels(), rotation = 270, horizontalalignment = "left")
ax.set_xlabel("core age")
ax.set_title("LINSIGHT score by FANTOM architecture age")
ax.legend(bbox_to_anchor=(1.5, 1.0))
plt.savefig("%sfantom_linsight_architecture_coremrca.pdf" %(RE), bbox_inches = "tight")


# In[25]:


from matplotlib import gridspec

order = ["simple", "complexenh"]
order1 = ["simple", "complex"]

fig = plt.figure(figsize = (12, 8))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 3])
ax0 = plt.subplot(gs[0])

sns.barplot(x = "arch", y = "linsight_score", data = sample,
#            showfliers = False, notch = True,
            order = order1,
           palette = pal,
           ax = ax0)

ax0.set_xlabel("")
ax0.set_xticklabels("")
ax0.set_ylabel("LINSIGHT score")
ax0.set_ylim(0,0.67)
sns.set("poster")
sample2 = base_df.sample(frac =0.03)
ax2 = plt.subplot(gs[1])
sns.barplot(x = "taxon2", y = "linsight_score", data = sample2.sort_values(by="mrca_2"),
              hue = "arch",
              palette = pal,
              hue_order = order,
              #join = False,
              #dodge = 0.25,
             ax = ax2)

ax2.set_xticklabels(ax2.get_xticklabels(),rotation = 270, horizontalalignment = "left")
ax2.set_xlabel("architecture age")

ax2.legend().remove()
ax2.set_xlabel("")
ax2.set_xticklabels("")
#ax2.set_yticklabels("")
ax2.set_ylabel("")
ax2.set_ylim(0,0.67)
ax2.set_xlabel("")
sns.set("poster")
plt.savefig("%sFig3c-JOINT_fantom_linsight_architecture_mrca_enh.pdf" %(RE), bbox_inches = "tight")


# In[112]:


old.head()


# In[121]:


order = ["simple", "complexenh"]
fig, ax = plt.subplots(figsize= (8,8))
sns.boxplot( x = "arch", y = "linsight_score", data = sample2.loc[sample2.mrca_2> 0.175],
            notch = True, order = order, palette = pal)
old = sample2.loc[sample2.mrca_2> 0.175]
m, mp = stats.mannwhitneyu(old.loc[old.arch.str.contains("complex"), "linsight_score"],
                          old.loc[old.arch.str.contains("simple"), "linsight_score"])
print(m,mp)


# In[114]:


old.groupby("arch")["linsight_score"].mean()


# In[29]:


sample2.head()


# In[52]:


kw_list = []
for i in sample2.mrca_2.unique():
    mrca_list = sample2.loc[sample2.mrca_2 == i, "linsight_score"].to_list()
    kw_list.append(mrca_list)
from scipy.stats import mstats

args=[l for l in kw_list]
stats.mstats.kruskalwallis(*args)


# In[ ]:


core_age = df.groupby(["enh_id"])["mrca"].max().reset_index()
core_age = pandas.merge(core_age, syn_gen_bkgd) # add in taxon2
core_age = core_age[["enh_id", "mrca_2"]]
core_age.columns = ["enh_id", "mrca_2_core"]
core_age.sort_values(by="mrca_2_core").head()


# # Two-way ANOVA does not work because
# Independent observations - yes
# Residue distribution is normal - no (Jarque-bera)
# homogeneity of variance - no (Omnibus)

# In[55]:


sample2.head()


# In[108]:


import statsmodels.api as sm
from statsmodels.formula.api import ols
import statsmodels.stats.multicomp
sample3 = base_df.sample(frac =0.01)
# Fits the model with the interaction term
# This will also automatically include the main effects for each factor
model = ols('linsight_score ~ C(mrca_2)*C(arch)', sample3).fit()

# Seeing if the overall model is significant
print(f"Overall model F({model.df_model: .0f},{model.df_resid: .0f}) = {model.fvalue: .3f}, p = {model.f_pvalue: .4f}")

model.summary()


# In[110]:


subset = sample3[['mrca', 'arch',]].drop_duplicates()
tuples = [tuple(x) for x in subset.to_numpy()]
tuples

resid = model.resid
factor_groups = sample3.groupby(['mrca_2','arch'])

plt.figure(figsize=(6,6));
for values, group in factor_groups:
    i,j = values
    group_num = i*100  # for plotting purposes
    x = [group_num] * len(group)
    plt.scatter(x, resid[group.index],
            s=10)
plt.xlabel('Group');
plt.ylabel('Residuals');


# In[57]:


res = sm.stats.anova_lm(model, typ= 2)
res


# In[ ]:
