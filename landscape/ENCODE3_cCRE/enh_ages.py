import os, sys
import pandas as pd
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats


ENHPATH = "/dors/capra_lab/projects/enhancer_ages/encode/hepg2/data/dELS_combined/"
ENHFILE = "dELS_combined_enh_age_arch_summary_matrix.bed"

ENHF = os.path.join(ENHPATH, "breaks", ENHFILE)


SHUFPATH = os.path.join(ENHPATH, "shuffle/breaks")
SHUFFILE = "shuf-dELS_combined_enh_age_arch_summary_matrix.bed"

SHUFF = os.path.join(SHUFPATH, SHUFFILE)

#%%
cols = ["chr", "start", "end", "enh_id","id", "max_seg", "core_remodeling", "arch", "mrca"]
df = pd.read_csv(ENHF, sep = '\t', header = None, usecols =[0,1,2,3,4,5,6,7,8], names = cols)

df.head()

#%%
shufdf = pd.read_csv(SHUFF, sep = '\t', header = None, usecols =[0,1,2,3,4,5,6,7,8], names = cols)

shufdf.head()
#%%
print(df.max_seg.median()) # median number of segments ==1
print(shufdf.max_seg.median()) # median number of segments ==1

#%%
catdf = pd.concat([df, shufdf])
catdf[["max_seg", "core_remodeling"]] = catdf[["max_seg", "core_remodeling"]].astype(int)

#formatting
catdf.loc[catdf.max_seg <= 1, "core_remodeling"] = 0
catdf.loc[catdf.core_remodeling < 1 , "arch"] = "simple"
catdf["enh_len"] = catdf.end - catdf.start

#%% basic info
catdf.groupby(["id", "arch"])["enh_id"].count()

#dELS simple = 58.5%
15312/(15312+10880)

#5x shuffle = 56.1%
71103/(71103+55612)

catdf.groupby(["id", "arch"])["enh_len"].describe()
#%%
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

fig, ax = plt.subplots(figsize = (6,6))
x, y = "arch", "mrca"
hue = "id"
sns.boxplot(x = x, y=y, data = catdf, notch = True, hue = hue)
ax.legend(bbox_to_anchor = (1,1))
