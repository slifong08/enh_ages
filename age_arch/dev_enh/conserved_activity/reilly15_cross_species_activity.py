import glob
import os, sys
import pandas as pd
import seaborn as sns
from scipy import stats

colors = [ "amber", "faded green"]
palette = sns.xkcd_palette(colors)
sns.palplot(palette)
plt.rcParams.update({'font.size': 15})
sns.set_style("white")

#%%
RE = "/dors/capra_lab/projects/enhancer_ages/reilly15/results/pleiotropy/"
path = "/dors/capra_lab/projects/enhancer_ages/reilly15/data/multiintersect/"
f = "%strim-0.5_multiintersect_hu_count.bed" % outpath

df = pd.read_csv(f, sep= '\t', header = None)
df.columns = ["chr_enh","start_enh", "end_enh",	"enh_id", "fourth_col", "id",\
"core_remodeling","arch","seg_index","mrca","enh_len", "taxon", "mrca_2",\
"taxon2",	"mya",	"mya2",	"seg_den",	"datatype", "count_overlap"]

df.head()

df.loc[df.seg_index >=5, "core_remodeling"] = 0
df.loc[df.seg_index <5, "arch"] = "simple"
#%% remove human, primate specific sequences. Sequence must be as olds as euarchontaglires (common ancestor with mouse) in order to be evaluated here.
# removes 470 sequences

print(len(df))
#df = df.loc[df.mrca_2>0.126]
print(len(df))

#%% plot simple v. complex
order = ["simple", "complexenh"]
sns.boxplot(x= "arch", y = "count_overlap", data = df, order = order)

#%%


from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator


fig = plt.figure(figsize = (12, 8))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 3])
ax0 = plt.subplot(gs[0])

sns.barplot(x = "arch", y = "count_overlap", data = df,
            palette = palette, order = order,
            ax = ax0)
nsimple = len(df.loc[df.arch == "simple"])
ncomplex = len(df.loc[df.arch != "simple"])
labels = ["n = %s"%nsimple, "n = %s"%ncomplex, ]
ax0.set_xticklabels(labels, rotation = 90)
ax0.set(xlabel="", ylabel ="Number of Active Species", ylim=(0,2))
ax0.yaxis.set_major_locator(MultipleLocator(0.5))

sns.set("poster")

ax2 = plt.subplot(gs[1])
sns.barplot(x = "taxon2", y = "count_overlap", hue = "core_remodeling",
              data = df.sort_values(by = "mrca_2"),
                palette = palette,
            ax = ax2)

ax2.set_xticklabels(ax2.get_xticklabels(), rotation = 90)
sns.set("poster")
ax2.yaxis.set_major_locator(MultipleLocator(0.5))
ax2.set(ylabel="",  ylim=(0,2))
ax2.legend().remove()
plt.savefig("%sFigS3b-JOINT_barplot_reilly15_cross_species_overlap_x_mrca_2.pdf" % RE, bbox_inches = "tight" )

#%%
mwu, p = stats.mannwhitneyu(df.loc[df.arch == "simple", "count_overlap"],\
df.loc[df.arch != "simple", "count_overlap"])
print(mwu, p)


#%%
df.groupby("arch")["count_overlap"].mean()

#%%
RE
