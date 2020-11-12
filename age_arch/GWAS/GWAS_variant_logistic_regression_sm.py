
#%% In[1]:

from collections import Counter

import glob

import matplotlib.pyplot as plt

import numpy as np
from numpy import mean

import os, sys
import pandas
from scipy import stats
import seaborn as sns
from sklearn.model_selection import train_test_split
import statsmodels.api as sm
import statsmodels.formula.api as smf
import subprocess

colors = [ "amber", "faded green"]
palette = sns.xkcd_palette(colors)
sns.palplot(palette)
plt.rcParams.update({'font.size': 15})
sns.set_style("white")

import datetime
LAST_RUN = datetime.datetime.now()
TODAY = (datetime.date.today())
RE = "/dors/capra_lab/projects/enhancer_ages/gwas_catalog/results/for_publication/"

print("last run", LAST_RUN)


#%% pleiotropy data


multipath = "/dors/capra_lab/projects/enhancer_ages/fantom/data/multiintersect/trimmed/"
multifile = "%strimmed_all_fantom_enh_112_tissues_multiintersect_0.5_count.bed"%multipath


#%% Import species data


syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pandas.read_csv(syn_gen_bkgd_file, sep = '\t')
syn_gen_bkgd["mrca"] = syn_gen_bkgd["mrca"].round(3)
syn_gen_bkgd["mrca_2"] = syn_gen_bkgd["mrca_2"].round(3)
syn_gen_bkgd["mya"] = syn_gen_bkgd["taxon2"].apply(lambda x: x.split(" ")[-1])

syn_gen_bkgd.head()


#%% intersect with GWAS overlaps only


gwas_path = "/dors/capra_lab/projects/enhancer_ages/gwas_catalog/data/"
gwasF = "%sgwasCatalog_2019-09-24_hg19_unique_cleaned_LDEx_p5e-8.bed" % gwas_path
outF= "%senh_x_var/FANTOM_x_gwas19_ldex.bed" % gwas_path

cmd = "bedtools intersect -a %s -b %s -wao > %s" % (multifile, gwasF, outF)
os.system(cmd)

print(outF)

#%% open pleiotropy


multi = pandas.read_csv(outF, sep ='\t', header = None, usecols=[0,1,2,3,4,5,6,7,15])
multi.head()


multi.columns = ["chr_enh", "start_enh", "end_enh", "old_len", "core_remodeling",
"mrca_2", "datatype", "count_overlap", "gwas_overlap"] # rename columns


multi.head()

#%%

multi.groupby("gwas_overlap")["datatype"].count()
#%%
"""
gwas_overlap
0    30218
1     1248
"""
#%%
multi["gwas_overlap"].unique()

#%% LOGISTIC REGRESSION FUNCTION
def logistic_regression(test_model, df):

    model = smf.logit(test_model, data = df,  freq_weights= df.wts)
    results = model.fit()
    return results


#%% format dataframe


multi["constant"] = 1 # add intercept

multi[["old_len", "core_remodeling", "count_overlap" ]] = multi[["old_len", "core_remodeling", "count_overlap" ]].astype(int)
multi["mrca_2"] = multi["mrca_2"].astype(float)
multi["wts"] = 1
multi.loc[multi.gwas_overlap ==0, "wts"] = 0.0042
#%%

test_model = 'gwas_overlap ~ old_len + core_remodeling + mrca_2 + count_overlap + core_remodeling*count_overlap'
results = logistic_regression(test_model, multi, )
print(results.summary())



#%%
results.pred_table()


#%% UNDERSAMPLE
from sklearn.utils import shuffle


for i in np.arange(10):
    minority = multi.loc[multi.gwas_overlap == 1]
    num_minority = len(minority) # number of sampled negatives
    new_majority = multi.loc[multi.gwas_overlap != 1].sample(n = num_minority)
    undersampled_df = pd.concat([minority, new_majority]) # make the testdf


    undersample_results = logistic_regression(test_model, undersampled_df)
    print(undersample_results.summary())
    print(undersample_results.pred_table())
#%%
