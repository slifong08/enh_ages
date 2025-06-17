#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/bin/python
###
#   name    | mary lauren benton
#   created | 2017
#   updated | 2018.10.09
#           | 2018.10.11
#           | 2018.10.29
#           | 2019.02.01
#           | 2019.04.08
#
#   depends on:
#       BEDtools v2.23.0-20 via pybedtools
#       /dors/capra_lab/data/dna/[species]/[species]/[species]_trim.chrom.sizes
#
# /dors/capra_lab/users/bentonml/data/dna/[species]/[species]_blacklist_gap.bed

#update 2019-06-25 -
###


# In[2]:


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
from pybedtools.helpers import BEDToolsError, cleanup, get_tempdir, set_tempdir
from statsmodels.stats import multitest
from scipy import stats
from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator


import datetime
import seaborn as sns
import subprocess
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
greyish_green ="#82a67d"


last_run = datetime.datetime.now()
today = (datetime.date.today())
print("last run", datetime.datetime.now())

RE = "/dors/capra_lab/projects/enhancer_ages/gwas_catalog/results/"


# In[3]:


SPLIT_SHUF = 1
SPLIT_SHUF_ARCH = 1

TEST_FILENAME = "/dors/capra_lab/projects/enhancer_ages/gwas_catalog/data/gwasCatalog_2019-09-24_hg19_unique_cleaned_LDEx_p5e-8.bed"

ENH_PATH = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/download/h3k27ac_plus_h3k4me3_minus_peaks/Hsap_H3K27ac_plus_H3K4me3_minus_*/cleaned_Hsap_H3K27ac_plus_H3K4me3_minus_*/breaks/"
SHUF_PATH ="/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/download/h3k27ac_plus_h3k4me3_minus_peaks/Hsap_H3K27ac_plus_H3K4me3_minus_*/cleaned_Hsap_H3K27ac_plus_H3K4me3_minus_*/shuffle/breaks/"

#ENH_PATH = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/download/h3k27ac_plus_h3k4me3_minus_peaks/test/"
#SHUF_PATH = "%sshuf/" % ENH_PATH

num_threads = int(os.getenv('SLURM_CPUS_PER_TASK', 1))


# In[4]:


desc_df = pd.read_csv("/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/roadmap_hg19_sample_id_desc.csv", header = None)
desc_df.columns = ["sid", "desc"]
desc_df.head()


# # Evaluate Trimmed Roadmap enhancers

# In[5]:


enhFS = glob.glob("%s*summary_matrix.bed" % ENH_PATH)
len(enhFS)


# In[6]:


enh_dict = {}

enhFS = glob.glob("%s*summary_matrix.bed" % ENH_PATH)
for f in enhFS:
    fid = (f.split("/")[-1]).split(".")[0]

    sid = "E"+(fid.split("E")[1]).split("_")[0]

    enh_dict[sid] = f


# In[7]:


shuf_enh_dict = {}
list_sid = []

shufFS = glob.glob("%s*summary_matrix.bed" %SHUF_PATH)

for f in shufFS:
    fid = (f.split("/")[-1]).split(".")[0]
    sid = "E"+(fid.split("E")[1]).split("_")[0]

    shuf_enh_dict[sid] = f

    list_sid.append(sid) # append the sid to a list

shuf_enh_dict.keys()


# In[8]:


# how many shuffles don't have a matching enhancer age file?
for i in list_sid:
    if i not in enh_dict.keys():
        print("rerun", i)


# In[9]:


# how many shuffles don't have a matching enhancer age file?
no_shuffles_skip = []
for i in enh_dict.keys():
    if i not in list_sid:
        no_shuffles_skip.append(i)
        print("rerun shuffles", i)


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

    if exp_sum == 0: # prevent against zero overlaps, which skews CI interval calculations
        exp_sum += 1

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
    p_val = p_sum / nboot

    bs["p-value"] = p_val
    print("empirical p", p_val)

    return ci, bs , p_val


def get_results(f, gwas_file, shuf_dict, freq):

    # collection dictionaries

    f_dict = {}
    bs_dict = {}

    # enhancer file
    s = BedTool(f)

    # how many matched shuffles are there?
    iterations = len(shuf_dict)

    print('python {:s} {:s}'.format(' '.join(sys.argv), str(datetime.datetime.now())))

    # generate variables for annotating analyses
    sid = (f.split("/")[-1]).split("-")[1]
    arch = ((f.split("/")[-1]).split("-")[0])

    target = (gwas_file.split("/")[-1]).split(".")[0]
    key = "_".join([sid, arch])


    # intersect enhancer w/ GWAS variants

    obs_sum, obs_int = calculateObserved(s, BedTool(gwas_file))

    count = len(open(f).readlines(  )) # count the number of enhancers tested
    obs_freq = (obs_sum/count) # get % of enhancers overlapping GWAS variants


    # print the num of enhancers, num enhancer overlapping GWAS, and frequency of enhancers overlapping GWAS

    print(count, "obs_sum", obs_sum, "obs_freq", obs_freq)


    # create pool and intersect shuffled w/ GWAS overlaps in parallel

    pool = Pool(num_threads)
    partial_calcExp = partial(calculateExpected, BedTool(gwas_file), shuf_dict)
    exp_sum_list = pool.map(partial_calcExp, [i for i in range(iterations)])


    # to get freq of shuf overlapping GWAS, count total shufs
    # note - shuf_dict keys and exp_sum_list are indexed the same, thus, indexed values match.

    exp_count_list = []

    for i in range(iterations):
        f = shuf_dict[i]
        exp_count_list.append((len(open(f).readlines(  ))))

    # divide number of GWAS overlaps by total shufs per shuf dataset

    exp_freq_list =([i / j for i, j in zip(exp_sum_list, exp_count_list)]) # get expected frequencies


    print("median exp_freq", np.median(exp_freq_list)) # print the median exp. %


    # wait for results to finish before calculating p-value

    pool.close()
    pool.join()

    # remove iterations that throw bedtools exceptions

    final_exp_sum_list = [x for x in exp_sum_list if x >= 0]
    exceptions = exp_sum_list.count(-999)

    # calculate p value from obs v. exp freq GWAS overlaps
    # returns a list of fold-change observations, and list of shuffled fold-changes.

    if exceptions / iterations <= .1:

        if freq == 0:
            obs_list, fold_changes = calculateEmpiricalP(obs_sum, exp_sum_list, freq) # calculate expectation on count of overlaps

        elif freq == 1:
            obs_list, fold_changes = calculateEmpiricalP(obs_freq, exp_freq_list, freq) # calculate expectation on frequency of overlap


        sd = np.std(fold_changes) # calculate the fold-change sd

        obs_fold_change = obs_list[3]

        # calculate CI + empirical p-value from bootstrapped exp dist of fold changes
        # returns ci list [0.025, 0.975], bootstrap dataframe, and empirical p-val

        ci, bs_df, emp_pval = bootstrapCI(fold_changes, obs_fold_change)

        print("ci", ci)
        obs_list.extend(ci)

        # add more info to list
        info =[sd, arch, sid, target] # add info
        obs_list.extend(info)


        # make a dataframe of all this information
        header = ['Observed', 'Expected', 'permuted_StdDev',
                  'FoldChange_Med', 'p-value',  'date_time',"ci_975", "ci_025",
                  "fold_change_std", "arch", "sid", "target"]

        enh = pd.DataFrame(data = [obs_list], columns = header)
        enh["arch"] = arch
        enh["p-value"] = emp_pval # replace  obs v. expect p-value with empirical p-value from bootstrap dist.

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


# In[11]:


def split_files(sid, testdir):

    # split up enh by architecture
    os.chdir(testdir)
    fid = "%s-TARGET" % sid
    cmd = '''awk '{print >$7"-%s.bed"}' %s ''' % (fid, enhF) # split up the files by simple/complex
    os.system(cmd)

    TARGET_FS = glob.glob("*TARGET.bed")

    return TARGET_FS


# In[12]:


def plot_results(df, sid):
    simple = df.loc[df.arch.str.contains("simple"), ["arch",  "FoldChange_Med", "yerr",
                                         "log2_FCmed", "log2yerr", "Observed",
                                         "fold_change_std"]]
    simMed, simErr = simple.FoldChange_Med, simple.yerr

    complexenh= df.loc[df.arch.str.contains("complex"), ["arch",  "FoldChange_Med", "yerr",
                                                 "log2_FCmed", "log2yerr", "Observed",
                                                "fold_change_std"]]
    comMed, comErr = complexenh.FoldChange_Med, complexenh.yerr

    sids = simple.arch

    ind = np.arange(len(comMed))  # the x locations for the groups
    width = 0.2  # the width of the bars
    barWidth = 0.25


    fig, ax = plt.subplots(figsize = (6,6))
    # Set position of bar on X axis
    r1 = np.arange(len(comMed))
    r2 = [x + barWidth for x in r1]


    # Make the plot
    plt.bar(r1, simMed, color=amber, width=barWidth, edgecolor='white', label='simple', yerr =simErr)
    plt.bar(r2, comMed, color=faded_green, width=barWidth, edgecolor='white', label='complexenh', yerr = comErr)

    result, p = stats.ttest_ind_from_stats(mean1 = simple.FoldChange_Med.item(), std1 =simple.fold_change_std.item(), nobs1 = simple.Observed.item(),
                mean2 = complexenh.FoldChange_Med.item(), std2 = complexenh.fold_change_std.item(), nobs2 = complexenh.Observed.item(),
                                       equal_var = False)
    plt.xlabel("%s, p = %s" % (sid,p))
    plt.ylabel("Fold-change")
    plt.xticks([r + barWidth for r in range(len(comMed))], sids, fontsize = 14)

    from matplotlib.ticker import MultipleLocator
    import matplotlib.ticker as ticker
    #ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(2**x))

    #ax.yaxis.set_major_formatter(ticks)
    for p in ax.patches:
        ax.annotate("%.2fx" % p.get_height(), (p.get_x() + p.get_width() / 2., p.get_height()-0.05),
                 ha='left', va='bottom', color='gray', xytext=(0, 10),
                 textcoords='offset points')

    # Create legend & Show graphic
    plt.legend(bbox_to_anchor = (1.45,1))

    ax.yaxis.set_major_locator(MultipleLocator(1))

    sns.set("poster")

    plt.savefig("%sfig4-ROADMAP_%s_matched_GWAS_2019_LDex_p5e-8_untrimmed.pdf"%(RE, sid))
    plt.show()


# In[13]:


# if relative simple needs to be assigned, reassign it according to the enhancer.

# split one shuffle file with 100 shuffles randomly into 100 smaller files

def reassign_relative_simple(shufF, enhF, sid):

    enh_count = len(open(enhF).readlines(  )) # get the count of enhancers

    shufP = "/".join(shufF.split("/")[:-1]) +"/"# get the shuf path

    relsimp = pd.read_csv(enhF, sep ='\t', header = None, usecols = [5])
    relative_simple = relsimp.median().iloc[0]

    relsimp_shuf = pd.read_csv(shufF, sep ='\t', header = None, usecols = [5])
    relative_simple_shuf = relsimp_shuf.median().iloc[0]

    # reassign shuffle simple/complex based on median number of enhancer age breaks

    if relative_simple !=relative_simple_shuf:
        print("reassign shuffled relative_simple from", relative_simple_shuf, relative_simple)

        if relative_simple > relative_simple_shuf:
            change_val_arg = '''$6<%d {$7="0"}1''' % relative_simple # turn add more simple enhancers to shuffles

        elif relative_simple < relative_simple_shuf:
            change_val_arg = '''$6>=%d {$7="1"}1'''% relative_simple # add more complex enhancers to shuffles

        temp = "%sshuf-temp.bed" % shufP
        cmd = '''awk 'BEGIN{FS=OFS="\t"} %s' %s > %s''' % (change_val_arg, shufF,temp) #& mv %s %s, temp, shufF)

        subprocess.call(cmd, shell = True)

        return temp

    else:
        print("relative simple is the same", relative_simple_shuf, relative_simple)
        return shufF

def split_shuffle(shufF, enhF, sid, testdir):

    # split up shuffles
    enh_count = len(open(enhF).readlines(  )) # get the count of enhancers

    for i in range(100):
        new_shuf = "%sshuf-%s-%d.bed" % (testdir, sid, i)
        cmd = "shuf -n %d %s > %s" % (enh_count, shufF, new_shuf)

        subprocess.call(cmd, shell = True)


# In[20]:


sid = "E044"
enhF = enh_dict[sid]
shufF = shuf_enh_dict[sid]

reassign_relative_simple(shufF, enhF, sid)


# In[14]:


CALCULATE_FOLD_CHANGE = 0
collection_dict = {}


# # SARAH - add already run files to collection_dict!

# In[15]:


FREQ = 1
for sid in enh_dict.keys():

    # output file + directory to temporarily hold files for analysis
    result_f = '/dors/capra_lab/projects/enhancer_ages/gwas_catalog/data/roadmap/figS4-ROADMAP_%s_GWAS_2019_LDex_p5e-8_untrimmed.txt' %sid
    testdir = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/download/h3k27ac_plus_h3k4me3_minus_peaks/Hsap_H3K27ac_plus_H3K4me3_minus_%s/cleaned_Hsap_H3K27ac_plus_H3K4me3_minus_%s/breaks/%s/" % (sid, sid, sid) # make a dir that you can delete.

    # get the corresponding files
    if (CALCULATE_FOLD_CHANGE == 1) & (sid not in collection_dict.keys()) & (sid not in no_shuffles_skip):
        enhF = enh_dict[sid] # get the enhancer file

        shufF = shuf_enh_dict[sid] # get the 100x shuf file

        # apply relative simple definition (from enhancer dataset) to shuffles
        reassigned_shuf = reassign_relative_simple(shufF, enhF, sid)


        # make the test directory
        cmd = "mkdir %s" %testdir
        os.system(cmd)

        os.chdir(testdir)
        print(sid, testdir)

        ### PROCESS THE MATCHED SHUFFLES ###

        # STEP 1 - split up shuffles by shuffle
        # split the shuffles.

        split_shuffle(reassigned_shuf, enhF, sid, testdir)


        # STEP 2 - split up shuffles by architecture

        shuffles= glob.glob("%s*shuf*.bed" % (testdir))

        for f in shuffles:

            fid = (f.split("/")[-1]).split(".")[0]

            os.chdir(testdir)
            cmd = '''awk '{print $1 "\t" $2 "\t" $3 >$7"-%s.bed"}' %s ''' % (fid, f) # split up the files by simple/complex
            os.system(cmd)


            cmd = "rm %s" % f # get rid of the input file.
            #os.system(cmd)

        ### SPLIT TRIMMED ENHANCER FILE on archtiecture###

        TARGET_FS = split_files(sid, testdir)

        simpleFS = glob.glob("%s0-*.bed" % (testdir)) # get all the simple shuffles

        complexFS = glob.glob("%s1-*.bed" % (testdir)) # get all the complex shuffles

        # get matched shuffled backgrounds

        simpleDict = dict((key, value) for key, value in enumerate(simpleFS)) # make a dictionary of simple files
        complexDict = dict((key, value) for key, value in enumerate(complexFS))  # make a dictionary of complex files

        results_dict ={}

        for sample in TARGET_FS:

            if "core_remodeling" not in sample:
                print(sample)

                if "0-" in sample:
                    print("simple")
                    shuf_arch_dict = simpleDict # get correct shuffle matched on architecture

                elif "1-" in sample:
                    print("complex")
                    shuf_arch_dict = complexDict # get correct shuffle matched on architecture

                resultsdf, expected_list = get_results(sample, TEST_FILENAME, shuf_arch_dict, FREQ) # bootstrap each of the files.
                results_dict[sample] = resultsdf

        df = pd.concat(results_dict.values())

        df["yerr"]=df["ci_025"] - df["ci_975"] # calculate yerr

        df["log2fold"] = np.log2(df.FoldChange_Med) # log2fold change of median enrichment

        val, df["fdr_p"] = multitest.fdrcorrection(df["p-value"], alpha = 0.05, method='indep') # FDR correction

        df.loc[df.sid.str.contains("0-"), "arch"] = "simple"
        df.loc[df.sid.str.contains("1-"), "arch"] = "complex"

        df.to_csv(result_f, sep = '\t', index = False, header = True) # write the file

        clean_up = "rm -r %s" % testdir

        os.system(clean_up)

        del results_dict


        collection_dict[sid] = df

        del df, sid

    elif os.path.exists(result_f) ==True:
        df= pd.read_csv(result_f, sep = '\t') # write the file

        collection_dict[sid] = df


        del df, sid


# In[17]:


CALCULATE_FOLD_CHANGE = 0
if CALCULATE_FOLD_CHANGE == 0:

    combined_dict = {}

    combined = glob.glob("/dors/capra_lab/projects/enhancer_ages/gwas_catalog/data/roadmap/figS4-ROADMAP_E*_GWAS_2019_LDex_p5e-8_untrimmed.txt")

    for f in combined:
        sid = (f.split("/")[-1]).split("_")[1]
        newdf = pd.read_csv(f, sep ='\t')
        combined_dict[sid]= newdf


    newdf = pd.concat(combined_dict.values())
    print(len(newdf))


# In[18]:


newdf.loc[newdf.arch == 0, "arch"] = "simple"
newdf.loc[newdf.arch == 1, "arch"] = "complex"
newdf.head()


# In[19]:



newdf = pd.merge(newdf, desc_df, how = "left", on = "sid")

newdf.head()


# In[20]:


# fdr<0.1
pal = sns.xkcd_palette(["amber", "faded green"])


# In[21]:


fdr_cutoff = 0.05
simple_count = 0
simple_sig = 0
complex_count = 0
complex_sig = 0

sig_list = []

for i in newdf.sid.unique():

    simplefc = newdf.loc[(newdf.sid == i) & (newdf.arch =="simple"), "FoldChange_Med"].iloc[0]

    simplep = newdf.loc[(newdf.sid == i) & (newdf.arch =="simple"), "p-value"].iloc[0]

    simpleci025 = newdf.loc[(newdf.sid == i) & (newdf.arch =="simple"), "ci_025"].iloc[0]
    simpleci975 = newdf.loc[(newdf.sid == i) & (newdf.arch =="simple"), "ci_975"].iloc[0]

    complexfc = newdf.loc[(newdf.sid == i) & (newdf.arch =="complex"), "FoldChange_Med"].iloc[0]

    complexp = newdf.loc[(newdf.sid == i) & (newdf.arch =="complex"), "p-value"].iloc[0]

    complexci025 = newdf.loc[(newdf.sid == i) & (newdf.arch !="simple"), "ci_025"].iloc[0]
    complexci975 = newdf.loc[(newdf.sid == i) & (newdf.arch !="simple"), "ci_975"].iloc[0]

    if simplefc>complexfc:
        simple_count+=1
        if complexci025 < simpleci975: # is the lower bound simple ci greater than the higher bound complex ci?
            simple_sig+=1
            sig_list.append(i)

    else:
        complex_count+=1
        if simpleci025 < complexci975: # is the lower bound complex ci greater than the higher bound simple ci?
            complex_sig+=1
            sig_list.append(i)

print(simple_count, complex_count, simple_sig, complex_sig)


# In[22]:


# all results
order = ["simple", "complex"]

fig, ax = plt.subplots(figsize = (6,6))
sns.boxplot(x = "arch", y = "FoldChange_Med", palette = pal,
              data = newdf, order = order, notch = True)
sns.swarmplot(x = "arch", y = "FoldChange_Med",
              data = newdf, order = order, linewidth = 1, palette = pal)

ax.legend(bbox_to_anchor = (1,1))
#ax.set_ylim(0,10)
#plt.savefig("%sfigS4.2-roadmap_matched_MRCA_fdr_all_box_matched_GWAS_2019_LDex_p5e-8.pdf"%(RE), bbox_inches = "tight")


# In[23]:


newdf.groupby("arch")["FoldChange_Med"].median()


# In[24]:


stats.mannwhitneyu(newdf.loc[(newdf.arch == "simple") , "FoldChange_Med"],
                  newdf.loc[(newdf.arch == "complex"), "FoldChange_Med"],)


# In[25]:


stats.mannwhitneyu(newdf.loc[(newdf.arch == "simple") & (newdf.fdr_p<fdr_cutoff), "FoldChange_Med"],
                  newdf.loc[(newdf.arch == "complex")& (newdf.fdr_p<fdr_cutoff), "FoldChange_Med"],)


# In[26]:


newdf.loc[newdf.fdr_p<fdr_cutoff].groupby("arch")["FoldChange_Med"].median()


# In[27]:


sig = newdf.loc[newdf.sid.isin(sig_list)]


# In[28]:


stats.mannwhitneyu(sig.loc[(sig.arch == "simple"), "FoldChange_Med"],
                  sig.loc[(sig.arch == "complex"), "FoldChange_Med"],)


# In[29]:


# fdr<0.1
fig, ax = plt.subplots(figsize = (6,6))
sns.pointplot(x = "arch", y = "FoldChange_Med",
              data = newdf.loc[newdf.sid.isin(sig_list)], hue = "desc", order = order,)
ax.legend(bbox_to_anchor = (1,1))
#ax.set_ylim(0,10)


# In[38]:



order = ["simple", "complex"]
fig, ax = plt.subplots(figsize = (6,6))
sns.set("poster")
test = newdf.loc[newdf.sid.isin(sig_list)]
sns.boxplot(x = "arch", y = "FoldChange_Med",
              data = test,
            notch = True,
            order = order, palette = pal)
sns.swarmplot(x = "arch", y = "FoldChange_Med",
              data = test,
              order = order, linewidth = 1, palette = pal)
sigN = len(test.sid.unique())

label = "significant enrichment\nn=%s/98 datasets" % sigN
ax.set(xlabel =label, ylabel = "GWAS catalog variant median fold enrichment")
#plt.savefig("%sfigS4.2-roadmap_matched_MRCA_fdr05_boxplot_matched_GWAS_2019_LDex_p5e-8_untrimmed.pdf"%(RE), bbox_inches = "tight")


# In[39]:


from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator
plotdf = newdf[["arch",  "FoldChange_Med", "yerr",
                "Observed", "sid",
                 "fold_change_std", "desc"]].sort_values(by = ["arch", "FoldChange_Med"], ascending = False)
simple = plotdf.loc[plotdf.arch== "simple"]
complexenh= plotdf.loc[plotdf.arch== "complex"]

simMeds = []
simErrs = []
comMeds = []
comErrs = []
sid_list = []

for i in plotdf.sid.unique():

    sid_list.append(desc_df.loc[desc_df.sid == i, "desc"].item())

    simple_test = simple.loc[simple.sid ==i]
    if len(simple_test)>0:

        simMed, simErr = simple_test.FoldChange_Med.item(), simple_test.yerr.item()
    else:
        simMed, simErr = 0,0
    simMeds.append(float(simMed))
    simErrs.append(simErr)


    complexenh_test = complexenh.loc[complexenh.sid ==i]
    if len(complexenh_test)>0:
        comMed, comErr = complexenh_test.FoldChange_Med.item(), complexenh_test.yerr.item()
    else:
        comMed, comErr = 0,0
    comMeds.append(comMed)
    comErrs.append(comErr)

sids = simple_test.arch
f, ax = plt.subplots(figsize = (8,40))
sns.set("poster")
ind = np.arange(len(comMeds))  # the x locations for the groups
width = 0.15  # the width of the bars
barWidth = 0.40

r1 = np.arange(len(comMeds))
r2 = [x + barWidth for x in r1]

# Set position of bar on X axis
# Make the plot
plt.barh(r1, simMeds, color=amber,
         height=barWidth,
         edgecolor='white', label='simple', xerr =simErrs)

plt.barh(r2, comMeds, color=faded_green,
         height=barWidth,
         edgecolor='white', label='complexenh', xerr = comErrs)
ax.barh(r1, simMeds, color=amber,
        height=barWidth,
        edgecolor='white', label='simple', xerr =simErrs)

ax.barh(r2, comMeds, color=faded_green,
        height=barWidth,
        edgecolor='white', label='complexenh', xerr = comErrs)

# Add xticks on the middle of the group bars

plt.xlabel("Fold-change")
plt.yticks([r + barWidth for r in range(len(comMeds))], sid_list)
ax.axvline(1, color = 'k', ls = "--")


ax.set_yticklabels(sid_list) #, rotation = 90)

#plt.savefig("%sfigS4.2-roadmap_matched_MRCA_all_fdr_matched_GWAS_2019_LDex_p5e-8_untrimmed.pdf"%(RE), bbox_inches = "tight")

plt.show()


# In[40]:


from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator
plotdf = newdf.loc[newdf.sid.isin(sig_list)]
plotdf = plotdf[["arch",  "FoldChange_Med", "yerr",
                "Observed", "sid",
                 "fold_change_std", "desc"]].sort_values(by = ["arch", "FoldChange_Med"], ascending = False)
simple = plotdf.loc[plotdf.arch== "simple"]
complexenh= plotdf.loc[plotdf.arch== "complex"]

simMeds = []
simErrs = []
comMeds = []
comErrs = []
sid_list = []

for i in plotdf.sid.unique():

    sid_list.append(desc_df.loc[desc_df.sid == i, "desc"].item())

    simple_test = simple.loc[simple.sid ==i]
    if len(simple_test)>0:

        simMed, simErr = simple_test.FoldChange_Med.item(), simple_test.yerr.item()
    else:
        simMed, simErr = 0,0
    simMeds.append(float(simMed))
    simErrs.append(simErr)


    complexenh_test = complexenh.loc[complexenh.sid ==i]
    if len(complexenh_test)>0:
        comMed, comErr = complexenh_test.FoldChange_Med.item(), complexenh_test.yerr.item()
    else:
        comMed, comErr = 0,0
    comMeds.append(comMed)
    comErrs.append(comErr)

sids = simple_test.arch
f, ax = plt.subplots(figsize = (8,40))
sns.set("poster")
ind = np.arange(len(comMeds))  # the x locations for the groups
width = 0.15  # the width of the bars
barWidth = 0.40

r1 = np.arange(len(comMeds))
r2 = [x + barWidth for x in r1]

# Set position of bar on X axis
# Make the plot
plt.barh(r1, simMeds, color=amber,
         height=barWidth,
         edgecolor='white', label='simple', xerr =simErrs)

plt.barh(r2, comMeds, color=faded_green,
         height=barWidth,
         edgecolor='white', label='complexenh', xerr = comErrs)
ax.barh(r1, simMeds, color=amber,
        height=barWidth,
        edgecolor='white', label='simple', xerr =simErrs)

ax.barh(r2, comMeds, color=faded_green,
        height=barWidth,
        edgecolor='white', label='complexenh', xerr = comErrs)

# Add xticks on the middle of the group bars

plt.xlabel("Fold-change")
plt.yticks([r + barWidth for r in range(len(comMeds))], sid_list)
ax.axvline(1, color = 'k', ls = "--")


ax.set_yticklabels(sid_list) #, rotation = 90)


plt.savefig("%sfigS4.2-roadmap_matched_MRCA_sig_median_matched_GWAS_2019_LDex_p5e-8_untrimmed.pdf"%(RE), bbox_inches = "tight")

plt.show()


# In[ ]:


len(plotdf.sid.unique())


# In[ ]:


df["dif"] = 0
for sid in df.sid2.unique():
    test = df.loc[df.sid2 == sid]
    break


# In[ ]:


test.head()


# In[ ]:
