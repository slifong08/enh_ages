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
#%% Files
path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/non-genic/"

summaryEnh ="%sno-exon_all_fantom_enh/breaks/no-exon_all_fantom_enh_ages_enh_age_arch_summary_matrix.bed" % path


summaryShuf = "%sno-exon_shuffle_fantom_age_arch_summary_matrix.bed" % path

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

#%%
shuffle = pd.read_csv(summaryShuf, sep = '\t' )

shuffle.head()

enh = pd.read_csv(summaryEnh, sep = '\t' )

#%%

shuffle = shuffle[["enh_id", "core_remodeling", "mrca_2", "taxon2"]].drop_duplicates()
enh = enh[["enh_id", "core_remodeling", "mrca_2"]].drop_duplicates()

shuffle

#%%


def fet_age(mrca, enh, shuffle, arch):

    # subset dataframes
    in_age_enh = enh.loc[enh.mrca_2 == mrca]
    in_age_shuf = shuffle.loc[shuffle.mrca_2 == mrca]


    # get counts
    in_arch = in_age_enh.loc[in_age_enh.core_remodeling==arch].shape[0]
    not_in_arch = in_age_enh.loc[in_age_enh.core_remodeling!=arch].shape[0]
    shuf_in_arch = in_age_shuf.loc[in_age_shuf.core_remodeling==arch].shape[0]
    shuf_not_in_arch = in_age_shuf.loc[in_age_shuf.core_remodeling!=arch].shape[0]

    # assign 2x2
    a = in_arch
    b = not_in_arch
    c = shuf_in_arch
    d = shuf_not_in_arch


    obs = [[a,b],[c,d]]
    print(obs)
    OR, P = stats.fisher_exact(obs)
    table = sm.stats.Table2x2(obs) # get confidence interval
    odds_ci = table.oddsratio_confint()
    newdf = pd.DataFrame({"mrca_2":[mrca], "a":[a], "b":[b], "c":[c], "d":[d],
                         "OR":[OR], "P":[P], "ci_lower" :[odds_ci[0]],
                        "ci_upper" :[odds_ci[1]], "core_remodeling_a":[arch]})

    print(mrca, obs, OR, P)
    return newdf

def fdr_correction(collection_dict):

    df = pd.concat(collection_dict.values())

    pvals = df["P"]

    df["rejected"], df["FDR_P"] = statsmodels.stats.multitest.fdrcorrection(pvals, alpha=0.05)
    return df
#%%


mrca_dict ={}
for mrca_2 in enh.mrca_2.unique():
    df = fet_age(mrca_2, enh, shuffle, 1)
    mrca_dict[mrca_2] = df

#%%

df = pd.concat(mrca_dict.values())
df = fdr_correction(mrca_dict)
df.sort_values(by = "mrca_2")
#%%
df["a"].sum()
df["b"].sum()
#%%
df.mrca_2 = df.mrca_2.round(3)
syn_gen_bkgd.mrca_2 = syn_gen_bkgd.mrca_2.round(3)
df = pd.merge(df, syn_gen_bkgd[["mrca_2", "taxon2"]], how = "left").drop_duplicates()

df["log2"] = np.log2(df["OR"])
df["yerr"] = df["ci_upper"] - df["ci_lower"]
df.sort_values(by = "mrca_2")


#%%

#%%
fig, ax = plt.subplots()
sns.set("poster")
sns.set_style("white")

x = "taxon2"
y ="log2"

data = df.loc[df.mrca_2>0].sort_values(by = "mrca_2")

splot = sns.barplot( x=x, y=y, data = data,
linewidth=2.5, facecolor=(1, 1, 1, 0), edgecolor=".2",
yerr =data["yerr"])

for n, p in enumerate(splot.patches):

    value = ORdf.iloc[n]["a"].astype(int)

    splot.annotate(value,
                   (p.get_x() + p.get_width() / 2.,0.2),
                   ha = 'center', va = 'baseline',
                   size=15,
                   rotation = 90,
                   color = "k",
                   xytext = (0, 1),
                   textcoords = 'offset points'
                   )
ax.set(ylabel= "Fold Change v. Bkgd\n(log2-scaled)",\
 title = "Complex enrichment per age", xlabel = "", ylim = (-2,1.2))

plt.axhline(0, color = "grey", linewidth = 2.5)

ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(round(2**x, 1)))

ax.yaxis.set_major_formatter(ticks)
ax.yaxis.set_major_locator(MultipleLocator(1))
ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
sns.set("poster")
sns.set_style("white")

plt.savefig("%sFigS6B_complex_odds_per_mrca.pdf" % RE, bbox_inches = 'tight')
