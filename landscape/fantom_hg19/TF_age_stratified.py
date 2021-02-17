import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns
import statsmodels
import statsmodels.api as sm

#%%
vals = {"path": "/dors/capra_lab/projects/enhancer_ages/fantom/data/tfbs/",
"syn_gen_bkgd_file": "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"}


vals["file"] = "%sall_fantom_enh_x_raw_tfbs_midpeak.bed" % vals["path"]


#%%


syn_gen_bkgd = pd.read_csv(vals["syn_gen_bkgd_file"], sep = '\t') # read the file
syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2", "mya", "mya2"]] # whittle down the df


#%%

# column names for dataframe
cols = ["chr_e", "start_e", "end_e", "enh_id",
"chr_s", "start_s", "end_s",
"seg_index", "core_remodeling", "core", "mrca",
"chr_t", "start_t", "enh_t", "tf", "reads?",
"midpeak", "raw_start_peak", "raw_end_peak",
"raw_len_peak", "mid_peak_len"]


df = pd.read_csv( vals["file"], sep = '\t', header = None) # open dataframe
df.columns = cols # name columns
df.mrca = df.mrca.round(3) # round these numbers
df = pd.merge(df, syn_gen_bkgd, how = "left", on = "mrca")


df.head()


#%% subset df by age
def tf_arch_odds(tf, df):

    collection = {}

    # get count of enhancer architectures with TF

    backbone = df.groupby(["mrca_2"])["enh_id"].count().reset_index()
    backbone.columns = ["mrca_2", "total_enh_count"]

    simpleTF = df.loc[(df.tf == tf) & (df.core_remodeling ==0)].groupby(["mrca_2"])["enh_id"].count().reset_index()
    simpleTF = pd.merge(backbone, simpleTF, how = "left").fillna(0)
    complexTF = df.loc[(df.tf == tf) & (df.core_remodeling ==1)].groupby(["mrca_2"])["enh_id"].count().reset_index()
    complexTF = pd.merge(backbone, complexTF, how = "left").fillna(0)
    print(simpleTF, complexTF)
    for mrca_2 in df.mrca_2.unique(): # stratify by age


        simple_age = simpleTF.loc[simpleTF.mrca_2 == mrca_2, "enh_id"].iloc[0]
        simple_not_age = simpleTF.loc[simpleTF.mrca_2 != mrca_2, "enh_id"].sum()

        complex_age = complexTF.loc[complexTF.mrca_2 == mrca_2, "enh_id"].iloc[0]
        complex_not_age = complexTF.loc[complexTF.mrca_2 != mrca_2, "enh_id"].sum()

        a, b, = simple_age, simple_not_age
        c, d = complex_age, complex_not_age

        obs = [[a,b], [c,d]] # does TF have preference for architecture and age?


        OR, P = stats.fisher_exact(obs)
        table = sm.stats.Table2x2(obs) # get confidence interval
        odds_ci = table.oddsratio_confint()
        newdf = pd.DataFrame({"mrca_2":[mrca_2], "a":[a], "b":[b], "c":[c], "d":[d],
                             "OR":[OR], "P":[P], "ci_lower" :[odds_ci[0]],
                            "ci_upper" :[odds_ci[1]], "TF":[tf],
                            "a_code": ["tf_simple_in_age"]})

        collection[mrca_2] = newdf

    returndf = pd.concat(collection.values())

    pvals = returndf["P"]

    returndf["rejected"], returndf["FDR_P"] = statsmodels.stats.multitest.fdrcorrection(pvals, alpha=0.05)

    return returndf


#%%
tf_results = {}

for tf in df.tf.unique():
    if tf != ".":
        results = tf_arch_odds(tf, df)

        tf_results[tf] = results
#%%


simple_age_enriched = pd.concat(tf_results.values())
simple_age_enriched["log2"] = np.log2(simple_age_enriched.OR)

#%%

table = pd.pivot(simple_age_enriched, values = "log2", index ="TF", columns = "mrca_2")
table = table.replace(-np.Inf, np.nan)
table = table.replace(np.Inf, np.nan)
#%%

table = table.sort_values(by = [0.00, 0.126],  ascending = False )


fig, ax = plt.subplots(figsize = (15,15))
sns.heatmap(table, center = 0, cmap="bwr")


""" interpretations -

This was a simple v. complex TFBS enrichment per age in FANTOM eRNA with ROADMAP ChIP-seq peaks.

Enrichment was measured as simple v. complex enhancers overlapping TF X in one age versus TF X in all other ages.

This answers the question - is a TF enriched for a specific age and architecture versus all other ages and architecture.

A limitation to this analysis is that TF-ChIP peaks are not matched to their cis-cell line.

Instead, they peaks that overlap any FANTOM eRNA locus.

ChIP-Peaks are 30bp long.

I don't know if I trust the results.
"""
