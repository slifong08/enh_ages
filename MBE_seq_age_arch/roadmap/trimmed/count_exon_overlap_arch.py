import glob
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

RE = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/results/for_publication/non-genic/trimmed/"
path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/download/h3k27ac_plus_h3k4me3_minus_peaks/trimmed/non-genic/"

noexon_fs = glob.glob("%sno-exon_trimmed*.bed" % path)
wexon_fs = glob.glob("%sexonOverlap_trimmed*.bed" % path)
print(len(noexon_fs))

def arch(fname, relative_simple):
    if "shuf" in fname:
        arch = pd.read_csv(fname, sep = '\t', header = None, usecols = [6])
    else:
        arch = pd.read_csv(fname, sep = '\t', header = None, usecols = [6])

    if relative_simple == 0:
        med = np.median(arch) # let non-exon file determine relative simple as 95% of enhancers do not overlap exons
    else:
        med = relative_simple


    simplen = len(arch.loc[arch[6]<=med])
    complexn = len(arch.loc[arch[6]>med])

    total = len(arch)
    print(med, simplen, complexn, total)
    return med, simplen, complexn, total
#%%
desc_path = "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19_roadmap_samples_enh_age/"
desc_df = pd.read_csv("%sroadmap_hg19_sample_id_desc.csv"% desc_path, header = None)
desc_df.columns = ["sid", "desc"]
desc_df.head()
#%%



sid_dict = {}
med_dict = {}


for f in noexon_fs:

    sid = (((f.split("/")[-1]).split("_")[1]).split("-")[-1]).split(".")[0]
    print(sid)

    if sid not in sid_dict:

        relative_simple = 0
        med, simplen, complexn, total = arch(f, relative_simple)

        newdf = pd.DataFrame({"median_breaks":[med],
        "simpleN": [simplen],
        "complexN": [complexn],
        "totalN": [total],
        "dataset":["no_ex"]
        })

        sid_dict[sid] = newdf # add
        med_dict[sid] = med


#%%
noexon_fs[0]
#%%

for f in wexon_fs:

    sid = (((f.split("/")[-1]).split("_")[1]).split(".")[0]).split("-")[-1]
    print(sid)

    if sid in sid_dict.keys():

        relative_simple = med_dict[sid]
        olddf = sid_dict[sid] # get the list

        med, simplen, complexn, total = arch(f, relative_simple)

        newdf = pd.DataFrame({"median_breaks":[med],
        "simpleN": [simplen],
        "complexN": [complexn],
        "totalN": [total],
        "dataset":["w_ex"]
        })

        df = pd.concat([olddf, newdf])

        sid_dict[sid] = df # add
#%%
print((sid_dict.keys()))



#%%
df = pd.concat(sid_dict.values())
df["frac_simple"] = df.simpleN.divide(df.totalN)
df["frac_complex"] = df.complexN.divide(df.totalN)
df.head()
#%%

df.sort_values(by = "frac_simple", ascending = False).head()



#%%
df.loc[df.dataset == "w_ex"].describe()
#%%

fig, ax = plt.subplots(figsize = (6,30))

sns.barplot(y= "sid2", x = "frac_ex",
data = df.sort_values(by = "frac_ex", ascending = False))
#ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
plt.savefig("%sexon_overlap_fraction_trimmed.pdf" %RE, bbox_inches = "tight")

#%%
fig, ax = plt.subplots(figsize = (6,30))

sns.barplot(y= "sid2", x = "ex_enh",
data = df.sort_values(by = "frac_ex", ascending = False))
#ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
plt.savefig("%sexon_overlap_count_trimmed.pdf" %RE, bbox_inches = "tight")


#%% match on the shuffled datasets you have available

limited = ['E038', 'E075', 'E022', 'E084', 'E067', 'E109', 'E117', 'E111', 'E129', 'E037', 'E043', 'E080', 'E087', 'E078', 'E085', 'E120', 'E079', 'E102', 'E101', 'E103', 'E021', 'E008', 'E099', 'E065', 'E092', 'E098', 'E097', 'E012', 'E090', 'E113', 'E089', 'E073', 'E061', 'E100', 'E114', 'E006', 'E104', 'E015', 'E095', 'E013', 'E121', 'E056', 'E014', 'E106', 'E125', 'E128', 'E026', 'E017', 'E105', 'E020', 'E049', 'E094', 'E068', 'E066', 'E096', 'E093', 'E007', 'E059', 'E112', 'E091', 'E011', 'E119', 'E005', 'E076', 'E003', 'E063', 'E004', 'E016']
lim = df.loc[df.sid.isin(limited)]

lim.describe()
#%%
lim.head()
#%%

fig, ax = plt.subplots(figsize = (6,10))

sns.barplot(y= "sid", x = "frac_ex",
data = lim.sort_values(by = "frac_ex", ascending = False))
#ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
plt.savefig("%sexon_overlap_fraction_limited.pdf" %RE, bbox_inches = "tight")

#%%
fig, ax = plt.subplots(figsize = (6,10))

sns.barplot(y= "sid", x = "ex_enh",
data = lim.sort_values(by = "frac_ex", ascending = False))
#ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
plt.savefig("%sexon_overlap_count_limited.pdf" %RE, bbox_inches = "tight")
