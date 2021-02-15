import os
import pandas as pd
from scipy import stats
import statsmodels
import statsmodels.api as sm
import subprocess

FANTOMPATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/hg38/HEPG2_FANTOM5_hg38/ages"
FANTOMFILE = "syn_breaks_HEPG2_FANTOM5_hg38_ages.bed"
FANTOM = os.path.join(FANTOMPATH, FANTOMFILE)

ENCODEPATH = "/dors/capra_lab/data/encode/encode3_hg38"
ENCODEFILE = "HepG2.bed"

ENCODE = os.path.join(ENCODEPATH, ENCODEFILE)

INTERSECTIONFILE = "HepG2_FANTOM_x_ENCODE_hg38.bed"
INTERSECTION = os.path.join(FANTOMPATH, INTERSECTIONFILE)


#%% Bed command

cmd = "bedtools intersect -a %s -b %s -wao > %s" % (FANTOM, ENCODE, INTERSECTION)

subprocess.call(cmd, shell = True)

print(cmd)


#%% dataframe


cols = ["chr", "start", "end",
"enh_id",
"chr_syn", "start_syn", "end_syn",
"seg_index", "core_remodeling", "core",
"mrca",
"chr_tf", "start_tf", "end_tf",
"tf_id", "reads", "tf_peak_centered_len",
"tf", "tf_peak_len", "cell_lines", "overlap"
]

df = pd.read_csv(INTERSECTION,
sep = '\t',
header = None,
names = cols)

df.info()


#%% remove chip-peaks with overlap less than 6bp

df = df.loc[df.overlap >5] # removes 1551 TF overlaps

df.info()


#%% add label to each enhancer.


df["arch"] = "simple"
df.loc[(df.core_remodeling ==1) & (df.core ==1), "arch"] = "complex_core"
df.loc[(df.core_remodeling ==1) & (df.core ==0), "arch"] = "complex_derived"


#%%


def prep_2x2(tf, arch, df):


    dfarch = df.loc[df.arch == arch]
    dfbkgd = df.loc[df.arch != arch]


    TF_in_arch = len(dfarch.loc[dfarch.tf == tf])
    TF_bkgd = len(dfbkgd.loc[dfbkgd.tf == tf])
    not_TF_in_arch = len(dfarch.loc[dfarch.tf != tf])
    not_TF_in_bkgd = len(dfbkgd.loc[dfbkgd.tf != tf])

    a, b, c, d = TF_in_arch, not_TF_in_arch, TF_bkgd, not_TF_in_bkgd

    obs = [[a,b], [c,d]]

    return obs


def quantify_2x2(obs, comparison_name):


    OR, P = stats.fisher_exact(obs)
    table = sm.stats.Table2x2(obs) # get confidence interval
    odds_ci = table.oddsratio_confint()
    newdf = pd.DataFrame({"comparison_name":comparison_name,
                          "a":obs[0][0], "b":obs[0][1],
                          "c":obs[1][0], "d":obs[1][1],
                          "OR":[OR], "P":[P],
                          "ci_lower" :[odds_ci[0]],
                          "ci_upper" :[odds_ci[1]],
                        })

    if P<0.05:
        print(comparison_name, obs, OR, P)

    return newdf


def fdr_correction(collection_dict):

    df = pd.concat(collection_dict.values())

    pvals = df["P"]

    df["reject_null"], df["FDR_P"] = statsmodels.stats.multitest.fdrcorrection(pvals, alpha=0.05)

    return df

#%%

collection_dict = {}

for tf in df.tf.unique():

    for arch in df.arch.unique():

        comparison_name = tf + "-" + arch
        obs = prep_2x2(tf, arch, df)

        results = quantify_2x2(obs, comparison_name)

        collection_dict[comparison_name] = results
#%%

results_df = fdr_correction(collection_dict)

results_df.loc[results_df.reject_null == True]
