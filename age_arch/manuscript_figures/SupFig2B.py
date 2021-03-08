import glob
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.ticker as ticker
import numpy as np
import os, sys
from scipy import stats
import seaborn as sns
import statsmodels
import statsmodels.api as sm

RE ="/dors/capra_lab/projects/enhancer_ages/fantom/results/for_publication/age_breaks/"

es_colors = [ "slate grey","greyish"]
es_pal = sns.xkcd_palette(es_colors)
sns.palplot(es_pal)

#%% Files
path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/"

enh = "%sFANTOM_enh_age_arch_full_matrix.tsv" % path
summaryEnh = "%sFANTOM_enh_age_arch_summary_matrix.tsv" % path

shuf = "%sSHUFFLED_FANTOM_enh_age_arch_full_matrix.tsv" % path
summaryShuf = "%sSHUFFLE_FANTOM_enh_age_arch_summary_matrix.tsv" % path

#%% other summary files

# age and taxon file
syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2", "mya", "mya2"]] # whittle down the df
syn_gen_bkgd

# tissue/cell line dataset descriptions
desc_file = "/dors/capra_lab/data/fantom/fantom5/facet_expressed_enhancers/sample_id_descriptions.txt"
desc_df= pd.read_csv(desc_file, sep = '\t', header = None)

#%% LOAD Files
enh = "%sFANTOM_enh_age_arch_full_matrix.tsv" % path
summaryEnh = "%sFANTOM_enh_age_arch_summary_matrix.tsv" % path

shuf = "%sSHUFFLED_FANTOM_enh_age_arch_full_matrix.tsv" % path
summaryShuf = "%sSHUFFLE_FANTOM_enh_age_arch_summary_matrix.tsv" % path

shuffle = pd.read_csv(shuf, sep = '\t')
final_merge = pd.read_csv(enh, sep = '\t')
#%%
enh_lens = final_merge.groupby(["enh_id", "core_remodeling","enh_len", "arch"])["mrca_2"].max().reset_index()
enh_lens.mrca_2 = enh_lens.mrca_2.round(3)

enh_lens = pd.merge(enh_lens, syn_gen_bkgd[[ "mrca_2","taxon2", "mya2"]], how = "left", on = "mrca_2")
enh_lens = enh_lens.drop_duplicates()
print(enh_lens.shape)

enh_lens["datatype"]="FANTOM"
enh_lens.head()
#%%

shuf_len = shuffle.groupby(["enh_id", "core_remodeling",\
 "enh_len", "shuf_id", "arch"])["mrca_2"].max().reset_index()

print(shuf_len.shape)

shuf_len["datatype"]="SHUFFLE"
shuf_len.mrca_2 = shuf_len.mrca_2.round(3)
shuf_len = pd.merge(shuf_len, syn_gen_bkgd[["mrca_2", "taxon2", "mya2"]],\
 how = "left", on = "mrca_2")

shuf_len = shuf_len.drop_duplicates()
print(shuf_len.shape)

#%%
x = len(enh_lens) # how many enhancers are there?
enh_len_freq = enh_lens.groupby(["mrca_2", "core_remodeling"])["enh_id"]\
.count().divide(x).round(3).reset_index()
enh_len_freq.columns = ["mrca_2", "core_remodeling", "enh_freq"]

y = len(shuf_len) # how many shuffled enhancers are there?
shuf_len_freq = shuf_len.groupby(["mrca_2", "core_remodeling"])["enh_id"]\
.count().divide(y).round(3).reset_index()
shuf_len_freq.columns = ["mrca_2", "core_remodeling", "shuf_freq"]

# calculat difference between exp and observed frequency
len_freq = pd.merge(enh_len_freq, shuf_len_freq, how = 'left' )
len_freq["dif"] = len_freq.enh_freq / len_freq.shuf_freq
len_freq

shuf_len = shuf_len.loc[~(shuf_len.arch.str.contains("complex")&(shuf_len.mrca_2 == 0.0))]
enh_lens = pd.merge(enh_lens, len_freq, how = "left")
shuf_len = pd.merge(shuf_len, len_freq, how = "left")
shuf_len.head()

#%% Supplemental Figure 1.2 - FANTOM/SHUFFLE enhancer age x length
len_concat = pd.concat([enh_lens, shuf_len])

len_concat["datatype2"] = len_concat.datatype + "-"+ len_concat.arch
len_concat["datatype2"].unique()


e_colors = [ "slate grey","greyish",]
e_pal = sns.xkcd_palette(e_colors)
hue_order = ['FANTOM', 'SHUFFLE']

order =["Simple", "Complexenh"]
sns.lmplot(y = "enh_len", x = "mya2", data = len_concat,\
            hue = "datatype",  palette = e_pal,
            hue_order = hue_order, scatter = False, x_estimator = np.median)
plt.savefig("%sfigS1.2b-LM_enh_len.pdf" % RE, bbox_inches = "tight")
