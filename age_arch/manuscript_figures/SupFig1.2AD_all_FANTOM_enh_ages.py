import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
RE ="/dors/capra_lab/projects/enhancer_ages/fantom/results/age_breaks/"
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
shuffle.mrca_2 = shuffle.mrca_2.round(3)

final_merge = pd.read_csv(enh, sep = '\t')
final_merge.mrca_2 = final_merge.mrca_2.round(3)

#%%
ANNOTATIONS = 0

shuffle_dist = shuffle.groupby(["enh_id", "shuf_id"])["mrca_2"].max().reset_index() # calculate age distribution shuffled enhancers by enh_id

shuffle_dist["mrca_2"] = shuffle_dist["mrca_2"].round(3) # round ages


shuffle_dist2 = shuffle_dist.groupby(["mrca_2", "shuf_id"])["enh_id"].count().reset_index()
shuffle_dist2.columns = ["mrca_2", "shuf_id", "mrca_count"]

totals = shuffle_dist2.groupby(["shuf_id"])["mrca_count"].sum().reset_index()
totals.columns = ["shuf_id", "id_totals"]

shuffle_dist2 = pd.merge(shuffle_dist2, totals, how = "left")
shuffle_dist2['freq'] = shuffle_dist2.mrca_count.divide(shuffle_dist2.id_totals)
shuffle_dist2["dataset"] = "shuffle"

#%%
plot_dist = final_merge.groupby("enh_id")["mrca_2"].max().reset_index() # calculate age distribution of FANTOM enhancers by enh_id

plot_dist["mrca_2"] = plot_dist["mrca_2"].round(3) # round ages

totals = plot_dist.groupby("mrca_2")["enh_id"].count().reset_index()
totals.columns = ["mrca_2", "mrca_count"]
totals['freq'] = totals.mrca_count.divide(totals.mrca_count.sum())

plot_dist = pd.merge(plot_dist, totals, how = "left", on = "mrca_2")
plot_dist["shuf_id"] = "FANTOM"
plot_dist["dataset"] = "FANTOM"
dist = pd.concat([shuffle_dist2, plot_dist]) # concat shuffle and fantom age distributions
dist.head()

#%%
dist = pd.merge(dist, syn_gen_bkgd, how = "left", on = "mrca_2")


# calculate MWU shuffle v. FANTOM age distributions (from each enhancer)

mwustat, mwupval = stats.mannwhitneyu(shuffle_dist2["mrca_2"].loc[shuffle_dist2["dataset"] =="shuffle"],
                                     plot_dist["mrca_2"].loc[plot_dist["dataset"] =="FANTOM"] )
dist = dist[["taxon2", "mrca_2", "shuf_id", "dataset", "freq"]].drop_duplicates()

dist.sort_values(by = "mrca_2").head()
#%%
# Fold change of fantom/shuffle age frequency
fold_change_dict = {}
for i in dist.shuf_id.unique():
    if i != "FANTOM":
        fold_change = dist.loc[(dist.shuf_id.str.contains(i))|\
        (dist.shuf_id == "FANTOM")].pivot(index="mrca_2",\
         columns='dataset')['freq'].reset_index()
        fold_change = pd.merge(fold_change, syn_gen_bkgd[["mrca_2", "taxon2"]], how = "left").drop_duplicates()
        fold_change["fold_change"] = np.log2(fold_change["FANTOM"].divide(fold_change["shuffle"]))
        fold_change_dict[i] = fold_change

fc = pd.concat(fold_change_dict.values())
#%% plot set up

colors = ["greyish", "slate grey",]
palette = sns.xkcd_palette(colors)
sns.set("poster")

# plot
mwustat, mwupval = stats.mannwhitneyu(shuffle_dist2["mrca_2"].loc[shuffle_dist2["dataset"] =="shuffle"],
                                     plot_dist["mrca_2"].loc[plot_dist["dataset"] =="FANTOM"] )
fig, (ax2, ax1) = plt.subplots(ncols = 2, figsize=(16,8))

# line plot
sns.barplot(x = "taxon2", y= "freq",  data = dist.sort_values(by="mrca_2"),\
              hue = "dataset", palette = palette, ax = ax2)

ax2.set_xticklabels(ax2.get_xticklabels(), rotation = 90, horizontalalignment="left")
ax2.set(ylabel = "Frequency", xlabel = "Taxon \n MWU = %s, pval = %s"\
 % (mwustat, mwupval), ylim = (0, 0.6), title = "FANTOM Enhancer Age Distribution",)

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

ax1.set(ylabel = "log2(Fold-Change)", xlabel = "Taxon \n MWU = %s, pval = %s"\
 % (mwustat, mwupval), ylim = (-2.7, 2), title = "FANTOM Enhancer Age Fold-change")

plt.savefig("%sfig1b_figS1.2a-FANTOM_WHOLE_ENHANCER_MRCA_DIST_TAXON2.pdf" % \
(RE), bbox_inches='tight')
#%%
enh_ages = final_merge.groupby(['enh_id'])["mrca_2"].max().reset_index()
enh_ages["data_type"] = "fantom"
shuf_ages = shuffle.groupby(['enh_id', 'shuf_id',])["mrca_2"].max().reset_index()
shuf_ages["data_type"] = "shuffle"

shuf_ages.head()
plot_ages = pd.concat([enh_ages, shuf_ages])

print(plot_ages.groupby("data_type")["mrca_2"].mean())
m, mp = stats.mannwhitneyu(plot_ages.loc[plot_ages.data_type == "shuffle", "mrca_2"],
                        plot_ages.loc[plot_ages.data_type == "fantom", "mrca_2"])
print(m, mp)

""" RESULTS

enhancer v shuffle ages
data_type MEANS
fantom     0.272919
shuffle    0.252228

mwu = 32896644783.0 p = 0.0
"""
