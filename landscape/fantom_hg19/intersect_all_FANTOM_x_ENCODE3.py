import os
import pandas as pd
from scipy import stats
import statsmodels
import statsmodels.api as sm
import subprocess

FANTOMPATH = "/dors/capra_lab/projects/enhancer_ages/fantom/data/all_fantom_enh/ages/"
FANTOMFILE = "syn_breaks_all_fantom_enh_ages.bed"

FANTOM = os.path.join(FANTOMPATH, FANTOMFILE)

ENCODEPATH = "/dors/capra_lab/data/encode/encode3_hg38/TF/"
ENCODEFILE = "trimmed_encRegTfbsClusteredWithCells.liftOver.to.hg19.bed"

ENCODE = os.path.join(ENCODEPATH, ENCODEFILE)

INTERSECTIONFILE = "All_FANTOM_x_ENCODE_hg19.bed"
INTERSECTION = os.path.join(FANTOMPATH, INTERSECTIONFILE)

RE = "/dors/capra_lab/projects/enhancer_ages/fantom/results/encode3/"

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
"tf_id", "peak_len", "reads",
"tf",  "overlap"
]

df = pd.read_csv(INTERSECTION,
sep = '\t',
header = None)

df.head()
df.columns = cols

df.info()


#%% remove chip-peaks with overlap less than 6bp

df = df.loc[df.overlap >5] # removes 34802 TF overlaps

df.info()

df.head()


#%% add label to each enhancer.


df["arch"] = "simple"
df.loc[(df.core_remodeling ==1) & (df.core ==1), "arch"] = "complex_core"
df.loc[(df.core_remodeling ==1) & (df.core ==0), "arch"] = "complex_derived"

df["syn_id"] = df.chr_syn + ":" + df.start_syn.map(str) + "-" + df.end_syn.map(str)
df["enh_len"] = df.end - df.start
df["syn_len"] = df.end_syn - df.start_syn

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

    if obs[0][0] > 100:

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
            print(comparison_name, obs, round(OR, 2), round(P, 4))

        return newdf


def fdr_correction(collection_dict):

    df = pd.concat(collection_dict.values())

    pvals = df["P"]

    df["reject_null"], df["FDR_P"] = statsmodels.stats.multitest.fdrcorrection(pvals, alpha=0.2)

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
#print(len(results_df))
results_df

results_df.loc[results_df.reject_null == True]
#%%


simple_collection_dict = {}

for tf in df.tf.unique():

    arch = "simple"

    comparison_name = tf + "-" + arch
    obs = prep_2x2(tf, arch, df)

    results = quantify_2x2(obs, comparison_name)

    simple_collection_dict[comparison_name] = results
#%%

results_df = fdr_correction(simple_collection_dict)

results_df
#%%

results_df.loc[results_df.reject_null == True]


#%%

tf_density = df.groupby(["syn_id", "syn_len", "arch"])["tf"].count().reset_index()
tf_density.head()
tf_density["tf_density"] = tf_density["tf"].divide(tf_density.syn_len)
tf_density

x = "arch"
y = "tf_density"
data = tf_density
sns.barplot(x, y, data = data)

simple_tfden = tf_density.loc[tf_density.arch == "simple", "tf_density"]
core_tfden = tf_density.loc[tf_density.arch == "complex_core", "tf_density"]
derived_tfden = tf_density.loc[tf_density.arch == "complex_derived", "tf_density"]

print("simple v. core", stats.mannwhitneyu(simple_tfden, core_tfden))
print("core v. derived", stats.mannwhitneyu(derived_tfden, core_tfden))
