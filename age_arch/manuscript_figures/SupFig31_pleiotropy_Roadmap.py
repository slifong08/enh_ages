#!/usr/bin/env python
# coding: utf-8

# In[1]:


import glob
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator
import numpy as np
import os, sys
from scipy import stats
import seaborn as sns
order = ["simple", "complex"]
get_ipython().run_line_magic('matplotlib', 'inline')


# sns colors
arch_colors = ["amber", "dusty purple", "windows blue","greyish"]
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

RE = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/results/for_publication/pleiotropy/"

frac = 0.5


# In[2]:


# path to files.

path ="/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/non-genic/"

# get multiintersected files - each has been trimmed to dataset mean,
#exons excluded, multiintersected with 97 other ROADMAP enhancers.

multi_fs = glob.glob("%smultiintersect/trim-E*_multiintersect_%s_count.bed"% (path, frac))

# raw enhancers have been aged, and number of age segments counted.


age_fs = glob.glob("%sbreaks/no-exon_*_parallel_breaks_enh_age_arch_summary_matrix.bed" % path)


# other annotations about MRCA

syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"

syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t')

syn_gen_bkgd["mrca"]=syn_gen_bkgd["mrca"].round(3)

syn_gen_bkgd["mrca_2"]=syn_gen_bkgd["mrca_2"].round(3)

syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2"]]


# other annotations about datasets descriptions
desc_df = pd.read_csv("/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/roadmap_hg19_sample_id_desc.csv", header = None)

desc_df.columns = ["sid", "desc"]

desc_df.head()


# In[9]:


# make dictionaries of sid + multiintersect or age files

multi = {}
age = {}


for f in multi_fs:

    sid = "E" + ((f.split("/")[-1]).split("E")[1]).split("_")[0]

    multi[sid] = f


for f in age_fs:

    sid = "E" + ((f.split("/")[-1]).split("E")[1]).split("_")[0]

    age[sid] = f


# In[10]:


len(age.keys()), len(multi.keys())


# In[12]:


results_data = {}
results_stats ={}



for key, multi_f in multi.items():

    # open the multiintersect df

    multi = pd.read_csv(multi_f, sep = '\t', header = None, usecols = [3, 4, 5, 6])

    m_cols = ["enh_id", "old_len", "mean_len", "count_overlap"]

    multi.columns = m_cols


    # open age file

    age_f = age[key]

    agedf = pd.read_csv(age_f, sep = '\t', header = None, usecols = [3, 8, 9, 10, 12, 13]) # open the count_overlap df

    a_cols  = ["enh_id",  "seg_index", "mrca", "old_len",  "mrca_2", "taxon2"]

    agedf.columns = a_cols

    # filter any enhancers longer than 10kb


    # merge
    enh = pd.merge(multi,agedf, how = "left", on = ["enh_id", "old_len"])

    # clean up long enhancers, and assign relative simple

    enh = enh.loc[enh.old_len<10000]

    relative_simple = enh.seg_index.median()

    print("relative simple # of age segments <", relative_simple)

    enh["arch"] = "simple"
    enh.loc[enh.seg_index >= relative_simple, "arch"] = "complex"


    DESC = desc_df.loc[desc_df.sid == key, "desc"].iloc[0]
    enh["desc"] = DESC
    results_data[DESC] = enh


    t, pval = stats.mannwhitneyu(enh["count_overlap"].loc[enh["arch"].str.contains("complex")],                                       enh["count_overlap"].loc[enh["arch"].str.contains("simple")])
    print(DESC, t, pval)


    medians = enh.groupby(['arch'])["count_overlap"].median().reset_index()
    medians["measure"] = "median"

    means = enh.groupby(['arch'])["count_overlap"].mean().reset_index()
    means["measure"] = "mean"

    measurements = pd.concat([medians, means])
    measurements["t-stat"],  measurements["p"],  measurements["sid"], measurements["desc"] = [t, pval, key, DESC]
    results_stats[key] = measurements


# In[13]:


results_data.keys()


# In[14]:


res = pd.concat(results_stats.values())
res.head()


# In[32]:


sns.set("poster")
fig, ax = plt.subplots(figsize = (8,40))
sns.barplot(y="desc", x = "count_overlap", data = res.sort_values(by = "count_overlap"), hue = "arch",
           palette = palette, hue_order = order)
ax.get_legend().remove()
ax.set_xlabel("Sample Count")

#ax.set_yticklabels(ax.get_xticklabels(), rotation = 270)

#plt.savefig("%spleiotropy-Roadmap_trimmed_%s.pdf"%(RE, frac), bbox_inches = "tight")


# In[20]:


order


# In[30]:


sns.set("poster")
fig, ax = plt.subplots(figsize = (8,8))
sns.boxplot(x="arch", y = "count_overlap", data = res, order = order,
           palette = palette, hue_order = order)
sns.swarmplot(x="arch", y = "count_overlap", data = res, order = order,
           palette = palette, hue_order = order)
#ax.get_legend().remove()
ax.set_ylabel("Sample Count")

#ax.set_yticklabels(ax.get_xticklabels(), rotation = 270)

#plt.savefig("%spleiotropy-Roadmap_trimmed_%s.pdf"%(RE, frac), bbox_inches = "tight")


# In[34]:


# fdr<0.1
fig, ax = plt.subplots(figsize = (6,6))
sns.pointplot(x="arch", y = "count_overlap",
              data = res, hue = "desc", alpha = 0.5)
ax.legend(bbox_to_anchor = (1,1)).remove()
#ax.set_ylim(0,10)


# In[35]:


results2 = pd.concat(results_data.values())


# In[ ]:


results2.mrca =results2.mrca.round(3)
results2 = pd.merge(results2, syn_gen_bkgd, how = "left")
results2["taxon"] = results2["taxon2"].str.split( n = 4, expand = True)


# In[81]:


sns.set("poster")


order = ["simple", "complex"]
order1 = ["simple", "complex"]

fig = plt.figure(figsize = (12, 8))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 3])
ax0 = plt.subplot(gs[0])

sns.barplot(x = "arch", y = "count_overlap", data = results2, palette = palette, order = order,

            ax = ax0)

ax0.yaxis.set_major_locator(MultipleLocator(4))
sns.set("poster")

nsimple = len(results2.loc[results2.arch == "simple"])
ncomplex = len(results2.loc[results2.arch != "simple"])
ax0_xlabels = ["simple\nn=%d" % nsimple, "complex\nn=%d" % ncomplex]

ax2 = plt.subplot(gs[1])

sns.barplot(x = "taxon", y = "count_overlap", hue = "arch",
              data = test.sort_values(by = "mrca"),
                palette = palette,
            ax = ax2)

ax2.set_xticklabels(ax2.get_xticklabels(), rotation = 90)
ax2.set(xlabel = "", ylabel = "", title = ("98 Roadmap Datasets"))

ax2.legend().remove()
sns.set("poster")

ax2.yaxis.set_major_locator(MultipleLocator(4))

ax2lim = ax2.get_ylim()

ax0.set(xlabel = "", ylabel = "Number of Tissue/Cell Datasets", ylim = ax2lim)
ax0.set_xticklabels(ax0_xlabels, rotation = 90)
sns.set("poster")


plt.savefig("%sSupFig3.2_pleiotropy-Roadmap_trimmed_mean_noexon_98_datasets_%s_mrcas.pdf"%(RE, frac), bbox_inches = "tight")


# In[80]:


results2.groupby("arch")["count_overlap"].mean()


# In[78]:



for desc in results2.desc.unique():
    sns.set("poster")
    test = results2.loc[results2.desc == desc]

    order = ["simple", "complex"]
    order1 = ["simple", "complex"]

    fig = plt.figure(figsize = (12, 8))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 3])
    ax0 = plt.subplot(gs[0])

    sns.barplot(x = "arch", y = "count_overlap", data = test, palette = palette, order = order,

                ax = ax0)

    ax0.yaxis.set_major_locator(MultipleLocator(4))
    sns.set("poster")

    nsimple = len(test.loc[test.arch == "simple"])
    ncomplex = len(test.loc[test.arch != "simple"])
    ax0_xlabels = ["simple\nn=%d" % nsimple, "complex\nn=%d" % ncomplex]

    ax2 = plt.subplot(gs[1])

    sns.barplot(x = "taxon", y = "count_overlap", hue = "arch",
                  data = test.sort_values(by = "mrca"),
                    palette = palette,
                ax = ax2)

    ax2.set_xticklabels(ax2.get_xticklabels(), rotation = 90)
    ax2.set(xlabel = "", ylabel = "", title = ("%s"% desc))

    ax2.legend().remove()
    sns.set("poster")

    ax2.yaxis.set_major_locator(MultipleLocator(4))

    ax2lim = ax2.get_ylim()

    ax0.set(xlabel = "", ylabel = "Number of Tissue/Cell Datasets", ylim = ax2lim)
    ax0.set_xticklabels(ax0_xlabels, rotation = 90)
    sns.set("poster")


    plt.savefig("%sSupFig3.2_pleiotropy-Roadmap_trimmed_%s_%s_mrcas.pdf"%(RE, desc, frac), bbox_inches = "tight")


# In[82]:


RE


# In[ ]:


test.groupby(["arch", "mrca_2"])["enh_id"].count()


# In[ ]:


cmd = "cp /dors1/capra_lab/users/fongsl/enh_age/scripts/roadmap/sample_overlap/analysis_trim_ROADMAP_brain_blood_cellline_pleiotropy.ipynb /home/fongsl/enhancer_ages/roadmap/bin"
os.system(cmd)


# In[ ]:


def plot_reg(taxon2, df):

    test = df.loc[df.taxon2 == taxon2] # get age-specfic dataframe

    fig, ( ax2) = plt.subplots(figsize = (8,8))

    # architecture-specific dataframe
    simple = test.loc[test.arch == "simple"]

    complexenh = test.loc[test.arch != "simple"]

    x = "count_overlap"
    y = "old_len"


    # do linear regression for simple enhancers in age - pleiotropy x length
    # get coeffs of linear fit

    slope, intercept, r_value, p_value, std_err = stats.linregress(simple[x],simple[y])


    # plot simple enhancers pleiotropy x length


    # plot regplot w/ linear regression annotation
    sns.regplot(x=x, y=y, data = simple,
                line_kws={'label':"y={0:.1f}x+{1:.1f}".format(slope,intercept)},
                color="y", ax = ax2,
                x_estimator=np.median)


    # plot complex enhancers pleiotropy x length
    if taxon2 != "Homo sapiens (0)":

        slope, intercept, r_value, p_value, std_err = stats.linregress(complexenh[x],complexenh[y])


        sns.regplot(x=x, y=y, data = complexenh,
                line_kws={'label':"y={0:.1f}x+{1:.1f}".format(slope,intercept)},
                color="g", ax = ax2,
                x_estimator=np.median)

        ax2.set(title = "%s" % taxon2,
                xlim = (0,complexenh.count_overlap.max()),
                xlabel = 'context overlap',
                ylabel = "enhancer length (bp)")

        ax2.set_xticks(np.arange(0, complexenh.count_overlap.max(), step = 10))
        ax2.legend()
    plt.tight_layout()
    return fig


# In[ ]:


taxon2 = "Vertebrata (615)"
plot_reg(taxon2, results2)


# In[ ]:


for taxon2 in results2.taxon2.unique():

    fig = plot_reg(taxon2, results2)
    plt.savefig("%spleiotropy_x_length_arch-trimmed_mean_len_ROADMAP_%s.pdf" %(RE, taxon2))


# In[ ]:


results2.head()


# In[ ]:


#%% LENGTH X pleiotropy


results2["rounded_len"] = results2.old_len.round(-2) # round to the nearest 100.
results2.head()


# In[ ]:


#%%
def plot_reg_len(taxon2, df):

    test = df.loc[df.taxon2 == taxon2] # get age-specfic dataframe

    fig, (ax2) = plt.subplots(figsize = (8,8))

    # architecture-specific dataframe
    simple = test.loc[test.arch == "simple"]

    complexenh = test.loc[test.arch != "simple"]

    y = "count_overlap"
    x = "enh_len"


    # do linear regression for simple enhancers in age - pleiotropy x length
    # get coeffs of linear fit

    slope, intercept, r_value, p_value, std_err = stats.linregress(simple[x],simple[y])


    # plot simple enhancers pleiotropy x length


    # plot regplot w/ linear regression annotation
    sns.regplot(x=x, y=y, data = simple, x_bins = 20,
                line_kws={'label':"y={0:.3f}x+{1:.3f}".format(slope,intercept)},
                color="y", ax = ax2,
                x_estimator=np.mean)


    # plot complex enhancers pleiotropy x length
    if taxon2 != "Homo sapiens (0)":

        slope, intercept, r_value, p_value, std_err = stats.linregress(complexenh[x],complexenh[y])


        sns.regplot(x=x, y=y, data = complexenh, x_bins = 20,
                line_kws={'label':"y={0:.3f}x+{1:.3f}".format(slope,intercept)},
                color="g", ax = ax2,
                x_estimator=np.mean)

        ax2.set(title = "%s" % taxon2,
                #ylim = (0,complexenh.count_overlap.max()),
                ylabel = 'context overlap',
                xlabel = "enhancer length (bp)")

        #ax2.set_xticks(np.arange(0, complexenh.rounded_len.max(), step = 100), rotation = 90)
        ax2.legend()
    plt.tight_layout()
    return fig
#%%
taxon2 = "Vertebrata (615)"
plot_reg_len(taxon2, results2)
#%%
for taxon2 in results2.sort_values(by = "mrca_2").taxon2.unique():
    print(taxon2)
    fig = plot_reg_len(taxon2, results2)
    sid = taxon2.split(" ")[0]
    plt.savefig("%strimmed_length_x_pleiotropy_arch-%s.pdf" % (RE, sid), bbox_inches = "tight")


# In[ ]:
