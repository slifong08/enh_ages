
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
RE = "/dors/capra_lab/projects/enhancer_ages/fantom/results/for_publication/"

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
