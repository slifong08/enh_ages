import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
RE ="/dors/capra_lab/projects/enhancer_ages/reilly15/results/age_arch/"


#%% Files


path = "/dors/capra_lab/projects/enhancer_ages/reilly15/data/breaks/non-genic/"

#enh = "%sHsap_brain_enhancers_reilly15_enh_age_arch_full_matrix.tsv" % path
summaryEnh = "%sno-exon_Hsap_brain_enhancers_reilly15_gapexcluded_parallel_breaks_enh_age_arch_summary_matrix.bed" % path

#shuf = "%sHsap_brain_enhancers_reilly15_negs_age_arch_full_matrix.tsv" % path
summaryShuf = "%sno-exon_Hsap_brain_enhancers_reilly15_gapexcluded_3000_1set_negs_parallel_breaks_enh_age_arch_summary_matrix.bed" % path

#%% other summary files

# age and taxon file
syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2", "mya", "mya2"]] # whittle down the df


#%% LOAD Files

cols = ["enh_id", "shuf_id", "core_remodeling", "arch", "taxon2", "mrca_2", ]
shuffle = pd.read_csv(summaryShuf, sep = '\t', header = None,
 usecols=[3,5,6,7, 11, 12])
shuffle.columns = cols

shuffle.mrca_2 = shuffle.mrca_2.round(3)
print(shuffle.shape)

#%%
final_merge = pd.read_csv(summaryEnh, sep = '\t', header = None,
 usecols=[3,4,5,6, 11, 12])
final_merge.columns = cols
final_merge.mrca_2 = final_merge.mrca_2.round(3)
print(final_merge.shape)

final_merge.head()
#%%
 shuffle["shuf_id"].unique()
#%%
ANNOTATIONS = 0
shuffle["shuf_id"] = shuffle["shuf_id"]
shuffle_dist = shuffle.groupby(["enh_id", "shuf_id"])["mrca_2"].max().reset_index() # calculate age distribution shuffled enhancers by enh_id

shuffle_dist["mrca_2"] = shuffle_dist["mrca_2"].round(3) # round ages
shuffle_dist.head()

#%%

shuffle_dist2 = shuffle_dist.groupby(["mrca_2", "shuf_id"])["enh_id"].count().reset_index()
shuffle_dist2.columns = ["mrca_2", "shuf_id", "mrca_count"]

stotals = shuffle_dist2.groupby(["shuf_id"])["mrca_count"].sum().reset_index()
stotals.columns = ["shuf_id", "id_totals"]


shuffle_dist2 = pd.merge(shuffle_dist2, stotals, how = "left")
shuffle_dist2['freq'] = shuffle_dist2.mrca_count.divide(shuffle_dist2.id_totals)
shuffle_dist2["dataset"] = "shuffle"
shuffle_dist2.head()


#%%


plot_dist = final_merge.groupby(["enh_id"])["mrca_2"].max().reset_index() # calculate age distribution shuffled enhancers by enh_id

plot_dist["mrca_2"] = plot_dist["mrca_2"].round(3) # round ages
plot_dist["shuf_id"] = "Reilly15"
plot_dist.head()

plot_dist2 = plot_dist.groupby(["mrca_2", "shuf_id"])["enh_id"].count().reset_index()
plot_dist2.columns = ["mrca_2", "shuf_id", "mrca_count"]
totals = plot_dist.groupby("shuf_id")["enh_id"].count().reset_index()
totals.columns = ["shuf_id", "id_totals"]
totals



plot_dist2 = pd.merge(plot_dist2, totals, how = "left")
plot_dist2['freq'] = plot_dist2.mrca_count.divide(plot_dist2.id_totals)

plot_dist2["dataset"] = "Reilly15"
plot_dist2


#%%

dist = pd.concat([shuffle_dist2, plot_dist2]) # concat shuffle and reilly age distributions
dist.head()

#%%


dist = pd.merge(dist, syn_gen_bkgd, how = "left", on = "mrca_2")

#%%
shuffle_dist2.head()
#%%
# calculate MWU shuffle v. reilly age distributions (from each enhancer)

mwustat, mwupval = stats.mannwhitneyu(shuffle_dist2["mrca_2"], plot_dist["mrca_2"])
print(mwustat, mwupval)
dist = dist[["taxon2", "mrca_2", "shuf_id", "dataset", "freq"]].drop_duplicates()

dist.sort_values(by = "mrca_2").head()
#%%
i
dist.loc[(dist.shuf_id == i)|(dist.shuf_id == "Reilly15")].pivot(index="mrca_2", columns='dataset')['freq'].reset_index()
#%%
# Fold change of reilly/shuffle age frequency
fold_change_dict = {}
for i in dist.shuf_id.unique():
    if i != "Reilly15":
        print(i)
        fold_change = dist.loc[(dist.shuf_id.str.contains(i))|\
        (dist.shuf_id == "Reilly15")].pivot(index="mrca_2",\
         columns='dataset')['freq'].reset_index()
        fold_change = pd.merge(fold_change, syn_gen_bkgd[["mrca_2", "taxon2"]], how = "left").drop_duplicates()
        fold_change["fold_change"] = np.log2(fold_change["Reilly15"].divide(fold_change["shuffle"]))
        fold_change_dict[i] = fold_change

fc = pd.concat(fold_change_dict.values())


#%% MWU of ages in enhancer and shuffle datasets


enh_ages = final_merge.groupby(['enh_id'])["mrca_2"].max().reset_index()
enh_ages["data_type"] = "Reilly15"
shuf_ages = shuffle.groupby(['enh_id', 'shuf_id',])["mrca_2"].max().reset_index()
shuf_ages["data_type"] = "shuffle"

shuf_ages.head()
plot_ages = pd.concat([enh_ages, shuf_ages])

print(plot_ages.groupby("data_type")["mrca_2"].mean())
m, mp = stats.mannwhitneyu(plot_ages.loc[plot_ages.data_type == "shuffle", "mrca_2"],
                        plot_ages.loc[plot_ages.data_type == "Reilly15", "mrca_2"])
print(m, mp)


#%%
""" RESULTS

enhancer v shuffle ages
data_type MEANS
fantom     0.436485
shuffle    0.392575

mwu = 851439314.0 p = 0.0
"""
#%% plot set up

colors = ["greyish", "slate grey",]
palette = sns.xkcd_palette(colors)
sns.set("poster")

# plot

fig, (ax2, ax1) = plt.subplots(ncols = 2, figsize=(16,8))

# line plot
sns.barplot(x = "taxon2", y= "freq",  data = dist.sort_values(by="mrca_2"),\
              hue = "dataset", palette = palette, ax = ax2)

ax2.set_xticklabels(ax2.get_xticklabels(), rotation = 90, horizontalalignment="left")
ax2.set(ylabel = "Frequency", xlabel = "Taxon \n MWU = %s, pval = %s"\
 % (m, round(mp, 3)), ylim = (0, 0.6), title = "Reilly15 Enhancer Age Distribution",)

ax2.legend(frameon = False)

# fold change plot
colors = ["slate grey"]
palette_fc = sns.xkcd_palette(colors)
sns.barplot(x = "taxon2", y= "fold_change",  data = fc.sort_values(by="mrca_2"),\
              palette = palette_fc, ax = ax1)
if ANNOTATIONS == 1:
    for p in ax1.patches:
        x = p.get_x().round(2)
        y = (p.get_height().round(2))
        ax1.annotate('{:.1f}'.format(p.get_height()), (x, y), color = "black", fontsize = 20)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation = 90, horizontalalignment="left")

ax1.set(ylabel = "log2(Fold-Change)",
 ylim = (-9, 1.2),
 title = "Reilly15 Enhancer Age Fold-change")

plt.savefig("%sfig1b_figS1.2a-Reilly15_WHOLE_ENHANCER_MRCA_DIST_TAXON2.pdf" % \
(RE), bbox_inches='tight')
