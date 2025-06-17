#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/bin/python
###
#   name    | sarah fong
#   created | 2020-07-09
#
#   depends on:
#       BEDtools v2.23.0-20 via pybedtools
#
#       previously generated shuffle that is aged, length matched, and architecture matched
###


# In[1]:


import os
import sys, traceback
import argparse
import datetime
import glob
import numpy as np
import pandas as pd
from functools import partial
import matplotlib.pyplot as plt
from multiprocessing import Pool
from pybedtools import BedTool
from statsmodels.stats import multitest
from scipy import stats
import subprocess

import datetime
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')

faded_green = "#7bb274"
slate_grey = "#59656d"
amber = "#feb308"
dusty_lavender = "#ac86a8"
dusty_purple = "#825f87"
windows_blue = "#3778bf"
greyish = "#a8a495"
greyish_blue = "#5e819d"
greyish_purple = "#887191"
dull_yellow = "#eedc5b"
greyish_green = "#82a67d"

last_run = datetime.datetime.now()
today = (datetime.date.today())
print("last run", datetime.datetime.now())

RE = "/dors/capra_lab/projects/enhancer_ages/gwas_catalog/results/matched_bkgd/"


# In[2]:


# inputs

RUN_BOOTSTRAP = 1
SPLIT_SHUF_ARCH =0

TEST_FILENAME = "/dors/capra_lab/projects/enhancer_ages/gwas_catalog/data/gwasCatalog_2019-09-24_hg19_unique_cleaned_LDEx_p5e-8.bed"

# more shuffles
shuf_more_path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/shuffle/breaks/"
shufFS_more = glob.glob("%s**/shuf*_summary_matrix.bed" %shuf_more_path, recursive = True)
print(len(shufFS_more))

# more shuffles
SHUF_PATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/non-genic/"
shufFS = glob.glob("%sno-exon_SHUFFLE_FANTOM_shuf-all_fantom_enh_112_tissues-*_age_breaks_summary_matrix.bed" % SHUF_PATH)

len(shufFS)


# In[3]:


fs = ['/dors/capra_lab/projects/enhancer_ages/fantom/data/architecture_coordinates/all_fantom/no-exon_FANTOM_simple_age_breaks.bed',
      '/dors/capra_lab/projects/enhancer_ages/fantom/data/architecture_coordinates/all_fantom/no-exon_FANTOM_complexenh_age_breaks.bed']


# In[4]:


cols = ["chr", "start", "end", "core_remodeling", "id", 'len']
for f in shufFS:
    df = pd.read_csv(f, sep ='\t', header = None, usecols = [0,1,2,4,15])
    df[16] = df[2] - df[1]
    df.columns = cols
    break
df.head()


# In[5]:


if SPLIT_SHUF_ARCH ==1:

    os.chdir(SHUF_PATH)

    # need to split up shuffles by architecture
    for f in shufFS:

        fid = (f.split("/")[-1]).split(".")[0] # file id

        cmd = '''awk '{print >$5"-%s.bed"}' %s ''' % (fid, f) # split up the files by simple/complex
        subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=subprocess.PIPE, shell=True)
    for f in shufFS_more:

        fid = (f.split("/")[-1]).split(".")[0] # file id

        cmd = '''awk '{print >$7"-%s.bed"}' %s ''' % (fid, f) # split up the files by simple/complex
        subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=subprocess.PIPE, shell=True)



# In[6]:


simpleFS = glob.glob("%s0-*.bed" % SHUF_PATH) # get all the simple shuffles

simpleDict = dict((key, value) for key, value in enumerate(simpleFS)) # make a dictionary of simple files

complexFS = glob.glob("%s1*.bed" % SHUF_PATH) # get all the complex shuffles

complexDict = dict((key, value) for key, value in enumerate(complexFS))  # make a dictionary of complex files


# In[7]:


len(complexDict.keys())


# In[10]:


###
#   functions
###


def loadConstants(species):  # note chrom.sizes not used in current implementation | 2018.10.29
    return {'hg19': ("/dors/capra_lab/users/fongsl/data/ensembl/hg19_blacklist_gap_ensemblexon.bed", "/dors/capra_lab/data/dna/human/hg19/hg19_trim.chrom.sizes"),
            'hg38': ("/dors/capra_lab/users/bentonml/data/dna/hg38/hg38_blacklist_gap.bed", "/dors/capra_lab/data/dna/human/hg38/hg38_trim.chrom.sizes"),
            'mm10': ("/dors/capra_lab/users/bentonml/data/dna/mm10/mm10_blacklist_gap.bed", "/dors/capra_lab/data/dna/mouse/mm10/mm10_trim.chrom.sizes")
            }[species]


def calculateObserved(annotation, test):


    obs_sum = 0

    obs_intersect = annotation.intersect(test, wo=True) # intersect obs regions w/ test file (e.g. enh x gwas snps)

    for line in obs_intersect:
        obs_sum += int(line[-1]) # sum snps in obs regions

    return obs_sum, obs_intersect


def calculateExpected(test, shuf_dict, idx):

    annotation = BedTool(shuf_dict[idx]) # get the shuffle file from the dictionary

    exp_sum = 0

    exp_intersect = annotation.intersect(test, wo=True) # intersect shuffled obs regions w/ test file

    for line in exp_intersect:
        exp_sum += int(line[-1]) # sum snps in shuffled exp regions

    return exp_sum


def calculateEmpiricalP(obs, exp_sum_list, freq):

    mu = np.median(exp_sum_list)
    sigma = np.std(exp_sum_list)
    dist_from_mu = [exp - mu for exp in exp_sum_list]
    p_sum = sum(1 for exp_dist in dist_from_mu if abs(exp_dist) >= abs(obs - mu))

    if freq !=1:
        val = 1

    elif freq ==1:
        val = 0.000001

    fold_change = (obs + val) / (mu + val) # fold change of the median

    fold_changes = list((obs + val) / (m + val) for m in exp_sum_list) # get all the fold changes to calculate ci

    p_val = (p_sum + val) / (len(exp_sum_list) + val)

    print("Frequency?", freq, "obs", obs, "mu", mu, "fold_change", fold_change, "raw fold change", obs/mu)

    info = [obs, mu, sigma,
    fold_change, p_val,
     str(datetime.datetime.now())]

    return info, fold_changes


def bootstrapCI(fold_changes_list, obs): # bootstrap C.I.s

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

    bs = pd.DataFrame(data = bs_means, index = np.arange(nboot), columns = ["bs_means"])

    bs["deltas"] = bs.bs_means - xbar

    bs = bs.sort_values(by = "deltas", ascending= False) # sort the list.

    low = bs.deltas.quantile(0.025)
    high = bs.deltas.quantile(0.975)

    ci = xbar - [high, low]

    p_sum = sum(1 for bs_mean in bs.deltas if abs(bs_mean) >= abs(obs-xbar)) # find how many deltas are greater than obs delta.
    p_val = p_sum / n

    print("empirical p", p_val)

    return ci, bs


def get_results(f, gwas_file, shuf_dict, freq):

    f_dict = {}
    bs_dict = {}

    s = BedTool(f) # the sample file

    iterations = len(shuf_dict)

    print('python {:s} {:s}'.format(' '.join(sys.argv), str(datetime.datetime.now())))

    sid = (f.split("/")[-1]).split("_")[2]
    arch = ((f.split("/")[-1]).split("_")[-1]).split(".")[0]

    target = (gwas_file.split("/")[-1]).split(".")[0]
    key = "_".join([sid, arch])


    # run initial intersection and save overlap counts
    obs_sum, obs_int = calculateObserved(s, BedTool(gwas_file))

    count = len(open(f).readlines(  ))
    obs_freq = (obs_sum/count)


    print(count, "obs_sum", obs_sum, obs_freq)


    # create pool and run simulations in parallel
    pool = Pool(num_threads)
    partial_calcExp = partial(calculateExpected, BedTool(gwas_file), shuf_dict)
    exp_sum_list = pool.map(partial_calcExp, [i for i in range(iterations)])


    exp_count_list = []

    for i in range(iterations):
        f = shuf_dict[i]
        exp_count_list.append((len(open(f).readlines(  ))))

    exp_freq_list =([i / j for i, j in zip(exp_sum_list, exp_count_list)]) # get the mean expected frequency


    print("exp_sum", exp_freq_list[0], np.median(exp_freq_list))


    # wait for results to finish before calculating p-value
    pool.close()
    pool.join()

    # remove iterations that throw bedtools exceptions
    final_exp_sum_list = [x for x in exp_sum_list if x >= 0]
    exceptions = exp_sum_list.count(-999)

    # calculate empirical p value
    if exceptions / iterations <= .1:

        if freq == 0:
            obs_list, fold_changes = calculateEmpiricalP(obs_sum, exp_sum_list, freq) # calculate expectation on count of overlaps

        elif freq == 1:
            obs_list, fold_changes = calculateEmpiricalP(obs_freq, exp_freq_list, freq) # calculate expectation on frequency of overlap


        sd = np.std(fold_changes) # calculate the fold-change sd

        obs_fold_change = obs_list[3]

        ci, bs_df = bootstrapCI(fold_changes, obs_fold_change)
        print("ci", ci)
        obs_list.extend(ci)

        info =[sd, arch, sid, target] # add info
        obs_list.extend(info)

        header = ['Observed', 'Expected', 'permuted_StdDev',
                  'FoldChange_Med', 'p-value',  'date_time',"ci_975", "ci_025",
                  "fold_change_std", "arch", "sid", "target"]

        enh = pd.DataFrame(data = [obs_list], columns = header)
        enh["arch"] = arch

        bs_df["arch"] = arch
        bs_df["key"] = key

        f_dict[key] = enh
        bs_dict[key] = bs_df

    else:
        cleanup()
        sys.exit(1)

    df = pd.concat(f_dict.values())
    bs_df = pd.concat(bs_dict.values())

    df["obs_freq"] = obs_freq # add the observed overlap frequency
    df["exp_freq"] = np.median(exp_freq_list)
    df["exp_sd"] = np.std(exp_freq_list)# add the mean expected overlap frequency

    return df, bs_df


# # Run

# In[11]:


num_threads = int(os.getenv('SLURM_CPUS_PER_TASK', 1))

FREQ = 1
if RUN_BOOTSTRAP == 1:
    samples_lumped = []

    for f in fs:
        if "0" not in str(f):
            samples_lumped.append(f) # create list of files to be bootstrapped.

    print("to be bootstrapped", samples_lumped)

    results_dict ={}
    bs_dict = {}

    for sample in samples_lumped:

        if "simple" in sample:
            print("simple")
            shuf_dict = simpleDict # get correct shuffle matched on architecture

        elif "complexenh" in sample:
            print("complex")
            shuf_dict = complexDict # get correct shuffle matched on architecture

        resultsdf, bs_df = get_results(sample, TEST_FILENAME, shuf_dict, FREQ) # bootstrap each of the files.

        #resultsdf["shuffle_fold_changes"] = ' '.join(map(str, fold_changes)) # turn the cold changes into a string

        results_dict[sample] = resultsdf
        bs_dict[sample] = bs_df

    df = pd.concat(results_dict.values())

    df["yerr"]=df["ci_025"] - df["ci_975"] # calculate yerr

    df["log2fold"] = np.log2(df.FoldChange_Med) # log2fold change of median enrichment

    val, df["fdr_p"] = multitest.fdrcorrection(df["p-value"], alpha = 0.05, method='indep') # FDR correction

    df.to_csv('/dors/capra_lab/projects/enhancer_ages/gwas_catalog/data/fig4a_FANTOM_GWAS_2019_LDex_p5e-8.txt',
             sep = '\t', index = False, header = True) # write the file

elif RUN_BOOTSTRAP == 0:
    df = pd.read_csv('/dors/capra_lab/projects/enhancer_ages/gwas_catalog/data/fig4a_FANTOM_GWAS_2019_LDex_p5e-8.txt',
             sep = '\t')
    df.pd.read_csv('/dors/capra_lab/projects/enhancer_ages/gwas_catalog/data/fig4a_FANTOM_GWAS_2019_LDex_p5e-8_bootstraps.txt',
             sep = '\t', index = False, header = True) # write the file
    df["yerr"]=df["ci_025"] - df["ci_975"]  # calculate yerr


# In[14]:


bs = pd.concat(bs_dict.values())
bs.head()


# In[22]:


# get the p-value for complex enhancers in a simple bootstrapped distribution
archs = ["complexenh", "simple"]

for arch in archs:
    print(arch)
    if arch == "complexenh":
        other_arch = "simple"

    else:
        other_arch = "complexenh"


    obs = df.loc[df.sid.str.contains(arch), "FoldChange_Med"].iloc[0] # get obs fold-change

    bootstrapped_means = bs.loc[bs.key.str.contains(other_arch), "bs_means"].to_list() # get bs distribution of other arch
    bootstrapped_dist =  bs.loc[bs.key.str.contains(other_arch), "deltas"].to_list() # get the deltas distribution
    n = len(bootstrapped_dist)
    x_bar = np.mean(bootstrapped_means) # mean of exp fold-change

    obs_dist = abs(obs - x_bar) # calculate the distance from the arch mean of the other arch distribution

    p_sum = sum(1 for bs_mean in bootstrapped_dist if abs(bs_mean) >= obs_dist) # how many bootstraps are farther away than the observed?

    p_val = p_sum / n
    print(arch, p_val)


# In[13]:


df.head()


# In[23]:


from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator

simple = df.loc[df.sid== "simple", ["sid",  "FoldChange_Med", "yerr", "permuted_StdDev",
                                    "ci_975", "ci_025",   "Observed",
                                     "fold_change_std"]]
simMed, simErr = simple.FoldChange_Med, simple.yerr

complexenh= df.loc[df.sid== "complexenh", ["sid",  "FoldChange_Med", "yerr", "permuted_StdDev",
                                      "ci_975", "ci_025", "Observed",
                                            "fold_change_std"]]
comMed, comErr = complexenh.FoldChange_Med, complexenh.yerr

sids = simple.sid


# In[29]:


ind = np.arange(len(comMed))  # the x locations for the groups
width = 0.2  # the width of the bars
barWidth = 0.25


fig, ax = plt.subplots(figsize = (6,6))
# Set position of bar on X axis
r1 = np.arange(len(comMed))
r2 = [x + barWidth for x in r1]


# Make the plot
sns.set("poster")
sns.set_style("white")
plt.bar(r1, simMed, color=amber, width=barWidth, edgecolor='white', label='simple', yerr =simple.yerr.iloc[0])
plt.bar(r2, comMed, color=faded_green, width=barWidth, edgecolor='white', label='complexenh', yerr =complexenh.yerr.iloc[0])

# Add xticks on the middle of the group bars

plt.ylabel("Fold-change")

#plt.xticks([r + barWidth for r in range(len(comMed))], sids, fontsize = 14)

from matplotlib.ticker import MultipleLocator
import matplotlib.ticker as ticker

for p in ax.patches:
    ax.annotate("%.2fx" % p.get_height(), (p.get_x() + p.get_width() / 2., p.get_height()-0.05),
             ha='left', va='bottom', color='gray', xytext=(0, 10),
             textcoords='offset points')

# Create legend & Show graphic
plt.legend(bbox_to_anchor = (1.45,1)).remove()

ax.yaxis.set_major_locator(MultipleLocator(0.1))
plt.ylim(0.9, 1.3)


plt.savefig("%sfig4-FANTOM_no-exon_matched_GWAS_2019_LDex_p5e-8_freq.pdf"%(RE))
plt.show()


# In[30]:


complexenh.yerr.iloc[0]


# In[31]:


# Welch's test to compare the means of two independent samples

# sample 1 = simple enhancers fold_change_median
# sample 2 = complex enhancers fold_change_median
# StdDev = the SD will be the 95% confidence intervals
# equal variance = False

result, p = stats.ttest_ind_from_stats(mean1 = simple.FoldChange_Med.iloc[0],
                                       std1 = simple.yerr.iloc[0],
                                       nobs1 = 1088,
                                      mean2 = complexenh.FoldChange_Med.iloc[0],
                                       std2 = complexenh.yerr.iloc[0],
                                       nobs2 = 1088,
                                       equal_var = False)
print(result, p)


# In[32]:


# Welch's test to compare the means of two independent samples

# sample 1 = simple enhancers fold_change_median
# sample 2 = complex enhancers fold_change_median
# StdDev = the SD of the fold changes over shuffle
# equal variance = False

result, p = stats.ttest_ind_from_stats(mean1 = simple.FoldChange_Med.iloc[0],
                                       std1 = simple.fold_change_std.iloc[0],
                                       nobs1 = 690,
                                      mean2 = complexenh.FoldChange_Med.iloc[0],
                                       std2 = complexenh.fold_change_std.iloc[0],
                                       nobs2 = 589,
                                       equal_var = False)
print(result, p)


# In[ ]:


complexenh.Observed.iloc[0]


# In[ ]:


complexenh


# In[ ]:


simple


# In[ ]:


palette = ["amber", "faded green"]
pal= sns.xkcd_palette(palette)
fig, ax = plt.subplots()
sns.barplot( x = "sid", y = "Observed", data = df, palette = pal,)
ax.axhline(df["exp_freq"].iloc[0], color = amber, linestyle = "--")
ax.axhline(df["exp_freq"].iloc[1], color = faded_green, linestyle = "--")
ax.legend().remove()
ax.set(ylabel = "% SNP overlap", xlabel ="")
plt.savefig("%sSNP_freq_FANTOM_no_exon_matched_bkgd.pdf" % RE, bbox_inches = 'tight')


# In[ ]:


df.to_csv('/dors/capra_lab/projects/enhancer_ages/gwas_catalog/data/fig4a_FANTOM_GWAS_2019_LDex_p5e-8.txt',
             sep = '\t', index = False, header = True) # write the file


# In[ ]:


(690+589)/30474


# In[ ]:
