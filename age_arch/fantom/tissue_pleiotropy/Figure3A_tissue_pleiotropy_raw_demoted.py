#!/usr/bin/env python
# coding: utf-8

#%% In[1]:


import glob
import pandas
import matplotlib.pyplot as plt
import numpy as np
import os, sys
from scipy import stats
import seaborn as sns
import statsmodels


colors = [ "amber", "faded green"]
palette = sns.xkcd_palette(colors)
sns.palplot(palette)
plt.rcParams.update({'font.size': 15})
sns.set_style("white")

import datetime
LAST_RUN = datetime.datetime.now()
TODAY = (datetime.date.today())
RE = "/dors/capra_lab/projects/enhancer_ages/fantom/results/for_publication/pleiotropy/"

print("last run", LAST_RUN)

path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/download/tfbs/"
samples = glob.glob("%s*_enh_tfbs_density.bed" % path)

multipath = "/dors/capra_lab/projects/enhancer_ages/fantom/data/multiintersect/"
multifile = "%sraw_FANTOM_enh_age_arch_summary_matrix_multiintersect_0.5_count.bed"\
%multipath


# #%% Import species data


syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pandas.read_csv(syn_gen_bkgd_file, sep = '\t')
syn_gen_bkgd["mrca"] = syn_gen_bkgd["mrca"].round(3)
syn_gen_bkgd["mrca_2"] = syn_gen_bkgd["mrca_2"].round(3)
syn_gen_bkgd["mya"] = syn_gen_bkgd["taxon2"].apply(lambda x: x.split(" ")[-1])

syn_gen_bkgd.head()


# # get enhancer file tissue/cell line descriptions

#%% In[3]:


desc_file = "/dors/capra_lab/data/fantom/fantom5/facet_expressed_enhancers/sample_id_descriptions.txt"
desc_df= pandas.read_csv(desc_file, sep = '\t', header = None)


# #%% Initial load and merge all FANTOM eRNA x TFBS datasets together
# ### Time consuming. Do this once, save concatenated file, and re-open with arch below

# ## eRNA overlap between FANTOM datasets
#
# ### 29444 unique enhancer IDs among 97 FANTOM datasets
#
# ### measure enhancer redundancy by counting the number of samples w/ same enh_id

#%% In[6]:


multi = pandas.read_csv(multifile, sep ='\t', header = None)
multi.head()


#%% In[7]:


multi.columns = ["chr_enh", "start_enh", "end_enh", "enh_id",
 "core_remodeling", "arch", "seg_index", "mrca",
"enh_len", "taxon", "mrca_2", "taxon2", "mya",
"mya2", "seg_den", "datatype", "count_overlap"] # rename columns



multi.head(10)


#%% In[15]:


df = multi.groupby([ "enh_id","arch" ])["count_overlap"].max().reset_index()

df = df.loc[df.arch != ""]
df.head()




#%% In[18]:


sns.distplot(df["count_overlap"].loc[df["arch"].str.contains("simple")], label = "simple")
sns.distplot(df["count_overlap"].loc[df["arch"].str.contains("complex")], label = "complex")
plt.legend()


#%% In[19]:


t, pval = stats.mannwhitneyu(df["count_overlap"].loc[df["arch"].str.contains("complex")],
                                   df["count_overlap"].loc[df["arch"].str.contains("simple")])
print(t, pval)

simple_med = df["count_overlap"].loc[df["arch"].str.contains("simple")].median()
complex_med = df["count_overlap"].loc[df["arch"].str.contains("complex")].median()
print( "simple, complex medians", simple_med, complex_med)



#%% In[22]:


multi.groupby("core_remodeling")["chr_enh"].count()


#%% In[23]:


counts = multi.groupby(["mrca_2", "core_remodeling", "taxon2"])["chr_enh"].count().reset_index()
counts["chr_enh"] = counts["chr_enh"].astype(str) # turn counts into a str
countLabel = counts.groupby(["taxon2", "mrca_2"])["chr_enh"].apply(','.join).reset_index()
countLabel["label"] = countLabel.taxon2 +"\n(" + countLabel.chr_enh +")"
labels = countLabel.sort_values(by = "mrca_2")["label"].tolist()
labels


#%% In[24]:



multi.groupby(["mrca_2","arch"])["count_overlap"].count()


#%% In[26]:



df.groupby(["arch"])["count_overlap"].mean()


#%% In[36]:


from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator

order = ["simple", "complexenh"]
order1 = ["simple", "complexenh"]

fig = plt.figure(figsize = (12, 8))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 3])
sns.set("poster")
ax0 = plt.subplot(gs[0])

sns.barplot(x = "arch", y = "count_overlap", data = df,
            palette = palette, order = order,
            #showfliers = False, notch = True,
            ax = ax0)
ax0.set(xlabel= "", ylabel= "Number of Tissue/Cell Datasets", ylim=(-1, 16))
ax0.set_xticklabels("")

ax0.yaxis.set_major_locator(MultipleLocator(4))


ax1 = plt.subplot(gs[1])
sns.barplot(x = "taxon2", y = "count_overlap", hue = "core_remodeling",
              data = multi.sort_values(by = "mrca"),
                palette = palette,
            ax = ax1)
ax1.set(xlabel= "", ylabel= "", ylim=(-1, 16))
ax1.set_xticklabels(ax1.get_xticklabels(), rotation = 90)


ax1.legend().remove()



ax1.yaxis.set_major_locator(MultipleLocator(4))

plt.savefig("%sFig3a-JOINT_barplot_fantom_sample_overlap_x_mrca_2.pdf" % RE,
bbox_inches = "tight" )


#%% In[42]:


from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator

order = ["simple", "complexenh"]
order1 = ["simple", "complexenh"]

fig = plt.figure(figsize = (12, 8))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 3])
ax0 = plt.subplot(gs[0])

sns.boxplot(x = "arch", y = "count_overlap", data = df,
            palette = palette, order = order,
            showfliers = False, notch = True,
            ax = ax0)
ax0.set_xlabel("")
ax0.set_xticklabels("")
ax0.set_ylabel("Number of Tissue/Cell Datasets")
ax0.set_ylim(-1, 28)
ax0.yaxis.set_major_locator(MultipleLocator(4))
sns.set("poster")

ax1 = plt.subplot(gs[1])
sns.boxplot(x = "taxon2", y = "count_overlap", hue = "core_remodeling",
              data = multi.sort_values(by = "mrca"),
                palette = palette,
                showfliers = False, notch = True,
            ax = ax1)

ax1.set_xticklabels(ax1.get_xticklabels(), rotation = 90)
ax1.set_xlabel("")

ax1.legend().remove()
sns.set("poster")
ax1.set_ylim(-1, 28)
ax1.yaxis.set_major_locator(MultipleLocator(4))
ax1.set_xticklabels("")
ax1.set_ylabel("")

sns.set("poster")
plt.savefig("%sFig3a-JOINT_boxplot_fantom_raw_length_sample_overlap_x_mrca_2.pdf" % RE, bbox_inches = "tight" )


#%% In[21]:


# simple
a = sns.jointplot(x = "mrca", y = "count_overlap", data = multi.loc[multi.core_remodeling ==0])
a.annotate(stats.pearsonr)


#%% In[22]:


# complex
b = sns.jointplot(x = "mrca", y = "count_overlap", data = multi.loc[multi.core_remodeling ==1])
b.annotate(stats.pearsonr)


# #Kruskal Wallis

#%% In[23]:


kw_list = []
simple_list = []
complex_list = []
for i in multi.mrca_2.unique():
    mrca_list = multi.loc[multi.mrca_2 == i, "count_overlap"].to_list()
    kw_list.append(mrca_list)

    mrca_simple = multi.loc[(multi.mrca_2 == i) & (multi.core_remodeling == 0),  "count_overlap"].to_list()
    simple_list.append(mrca_simple)

    mrca_complex = multi.loc[(multi.mrca_2 == i) & (multi.core_remodeling == 1), "count_overlap"].to_list()
    complex_list.append(mrca_complex)

from scipy.stats import mstats

args=[l for l in kw_list]
argsSimple=[l for l in simple_list]
argsComplex=[l for l in complex_list]
stats.mstats.kruskalwallis(*args)


#%% In[24]:


stats.mstats.kruskalwallis(*argsSimple)


#%% In[25]:


stats.mstats.kruskalwallis(*argsComplex)

#%%
multi.head()
#%% In[26]:

def plot_reg(taxon2, df):

    test = df.loc[df.taxon2 == taxon2] # get age-specfic dataframe

    fig, (ax2) = plt.subplots(figsize = (8,8))

    # architecture-specific dataframe
    simple = test.loc[test.arch == "simple"]

    complexenh = test.loc[test.arch != "simple"]

    x = "count_overlap"
    y = "enh_len"


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
#%%
taxon2 = "Vertebrata (615)"
plot_reg(taxon2, multi)
#%%
for taxon2 in multi.sort_values(by = "mrca_2").taxon2.unique():
    print(taxon2)


    fig = plot_reg(taxon2, multi)
    sid = taxon2.split(" ")[0]
    plt.savefig("%spleiotropy_x_length_arch-%s.pdf" % (RE, sid), bbox_inches = "tight")

#%% LENGTH X pleiotropy
multi.head()
multi["rounded_len"] = multi.enh_len.round(-2) # round to the nearest 100.
multi.head()
#%%
def plot_reg_len(taxon2, df):

    test = df.loc[df.taxon2 == taxon2] # get age-specfic dataframe

    fig, (ax2) = plt.subplots(figsize = (8,8))

    # architecture-specific dataframe
    simple = test.loc[test.arch == "simple"]

    complexenh = test.loc[test.arch != "simple"]

    y = "count_overlap"
    x = "rounded_len"


    # do linear regression for simple enhancers in age - pleiotropy x length
    # get coeffs of linear fit

    slope, intercept, r_value, p_value, std_err = stats.linregress(simple[x],simple[y])


    # plot simple enhancers pleiotropy x length


    # plot regplot w/ linear regression annotation
    sns.regplot(x=x, y=y, data = simple,
                line_kws={'label':"y={0:.3f}x+{1:.3f}".format(slope,intercept)},
                color="y", ax = ax2,
                x_estimator=np.median)


    # plot complex enhancers pleiotropy x length
    if taxon2 != "Homo sapiens (0)":

        slope, intercept, r_value, p_value, std_err = stats.linregress(complexenh[x],complexenh[y])


        sns.regplot(x=x, y=y, data = complexenh,
                line_kws={'label':"y={0:.3f}x+{1:.3f}".format(slope,intercept)},
                color="g", ax = ax2,
                x_estimator=np.median)

        ax2.set(title = "%s" % taxon2,
                ylim = (0,complexenh.count_overlap.max()),
                ylabel = 'context overlap',
                xlabel = "enhancer length (bp)")

        #ax2.set_xticks(np.arange(0, complexenh.rounded_len.max(), step = 100), rotation = 90)
        ax2.legend()
    plt.tight_layout()
    return fig
#%%
taxon2 = "Vertebrata (615)"
plot_reg_len(taxon2, multi)
#%%
for taxon2 in multi.sort_values(by = "mrca_2").taxon2.unique():
    print(taxon2)
    fig = plot_reg_len(taxon2, multi)
    sid = taxon2.split(" ")[0]
    plt.savefig("%sraw_length_x_pleiotropy_arch-%s.pdf" % (RE, sid), bbox_inches = "tight")
