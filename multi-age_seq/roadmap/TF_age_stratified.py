import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import seaborn as sns
from scipy import stats
import statsmodels
import statsmodels.api as sm
import subprocess

#%% paths, values, etc.


vals = {"chip_path": "/dors/capra_lab/data/encode/",
"enh_path": "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/download/h3k27ac_plus_h3k4me3_minus_peaks/age_breaks/",
"out_path": "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/tfbs/ucsc_tfbs/histone_x_chip/",
"desc_df": "/dors/capra_lab/projects/enhancer_ages/roadmap_encode/data/hg19/roadmap_hg19_sample_id_desc.csv",
"syn_gen_bkgd": "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"}

vals["enh"] = "%sno-exon_E118_age_breaks.bed" % vals["enh_path"] # H3K27ac+ H3K4me3- HepG2 enhancers, aged, stratified by syntenic blocks.
vals["chip"] = "%sHepG2_Hepatocellular_Carcinoma_encode_chipseq.bed" % vals["chip_path"] # Encode2 ChIP-seq data for HepG2

vals["outfile"] = "%sHepG2_histone_x_chip.bed" % vals["out_path"] # intersection .bed file


### sometimes the enh file has headers that need to be removed.
# also, there are duplicates in these files that need to be removed.


def remove_header(full_path):

    remove_header = "sed -i '/start_syn/d' %s" % full_path

    subprocess.call(remove_header, shell = True)
    print("removed header!")


def remove_duplicates(file, path):

    remove_duplicates = "sort -k1,1 -k2,2 -k3,3 %s | uniq > %st && mv t %s" \
    % (file, path, file)

    subprocess.call(remove_duplicates, shell = True)

    print("removed duplicates")


def format_dataframe(df):


    df = df.loc[df.chr_syn != "chrX"] # remove chromosome X
    df = df.loc[df.enh_len <10000] # remove enhancers longer than 10kb
    df = df.loc[df.len_overlap <6] # remove enhancers w/ little peak overlap
    return df


def apply_mrca_groups(df, syn_gen_bkgd_file):

    syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t')
    syn_gen_bkgd[["mrca", "mrca_2"]]=syn_gen_bkgd[["mrca", "mrca_2"]].round(3)
    syn_gen_bkgd = syn_gen_bkgd[["mrca", "mrca_2", "taxon2"]]

    df["mrca"] = df["mrca"].round(3)
    df = pd.merge(df, syn_gen_bkgd, how = "left", on = "mrca")

    return df



#%% Intersect enhancers and ChIP-seq data from matched cell_lines


def enh_x_chip_intersection(enh, chip, outfile):

    cmd = "bedtools intersect -a %s -b %s -wao > %s" % (enh, chip, outfile)

    subprocess.call(cmd, shell = True)


#%% Intersect enhancers and chip-seq


enh_x_chip_intersection(vals["enh"],vals["chip"], vals["outfile"])


#%% Read the file into a dataframe

df = pd.read_csv(vals["outfile"], sep = '\t', header = None).drop_duplicates()

# name the columns

cols = ["chr_syn", "start_syn", "end_syn", "enh_id",
        "chr_enh", "start_enh", "end_enh", "seg_index",
        "core_remodeling", "core", "mrca", "sample_id",
        "enh_len",
        "chr_tf", "midstart_tf", "midend_tf", "tf",
        "expt_score", "cell_line", "tf_id", "len_overlap"]
df.columns = cols

df.head()

#%% format the dataframe

df = format_dataframe(df)

df = apply_mrca_groups(df, vals["syn_gen_bkgd"])

#%% evaluate the complex enhancers


cmplx = df.loc[df.core_remodeling == 1]
cmplx.head()


#%% count tfs in cores, derived per age


counts = cmplx.groupby(["core", "mrca_2", "tf"])["enh_id"].count().reset_index()
counts.head()


#%% plot the cores


x = "enh_id"
y= "tf"
data = counts.loc[(counts.tf != ".") & (counts.core ==1)].sort_values(by = "enh_id")
hue = "mrca_2"

fig, ax = plt.subplots(figsize = (6,10))
sns.barplot(x = x, y = y,
 data = data, hue = hue,)

ax.set(xlabel = "count of cores", ylabel = "tf", title = "cores")


#%% plot the derived

x = "enh_id"
y= "tf"
data = counts.loc[(counts.tf != ".") & (counts.core ==0)]
hue = "mrca_2"

fig, ax = plt.subplots(figsize = (6,15))
sns.barplot(x = x, y = y,
 data = data, hue = hue,)

ax.set(xlabel = "count of derived", ylabel = "tf", title = "derived")

#%%

core_tf = "FOXA1"
derived_tf = "CEBPB"

#%%

def test_cooccurrence(core_tf, derived_tf, cmplx):
    all_enhs_w_tfs = cmplx.loc[(cmplx.tf ==core_tf) |( cmplx.tf == derived_tf), "enh_id"].unique()

    core_enh_ids = cmplx.loc[(cmplx.core ==1) &( cmplx.tf == core_tf), "enh_id"].unique()
    der_enh_ids = cmplx.loc[(cmplx.core ==0) &( cmplx.tf == derived_tf), "enh_id"].unique()

    together = len(np.intersect1d(core_enh_ids, der_enh_ids)) # find the core and derived TF in the same enhancer.
    not_together_core = len(core_enh_ids) - together
    not_together_der = len(der_enh_ids) - together
    not_together_config = len(all_enhs_w_tfs) - (together + not_together_core + not_together_der)

    a, b = together, not_together_der
    c, d = not_together_core, not_together_config

    obs = [[a,b], [c,d]]

    OR, P = stats.fisher_exact(obs)
    table = sm.stats.Table2x2(obs) # get confidence interval
    odds_ci = table.oddsratio_confint()
    if p<0.05:
        print(OR, P, obs)

    newdf = pd.DataFrame({"TF_core":[core_tf], "TF_der":[derived_tf],
        "a":[a], "b":[b], "c":[c], "d":[d],
        "OR":[OR], "P":[P], "ci_lower" :[odds_ci[0]],
        "ci_upper" :[odds_ci[1]]})

    return newdf

#%%

counts = cmplx.groupby(["core", "tf"])["enh_id"].count().reset_index()
core_tfs = counts.loc[(counts.core ==1) & (counts.enh_id>10) & (counts.tf != "."), "tf"]
der_tfs = counts.loc[(counts.core ==0) & (counts.enh_id>1)& (counts.tf != "."),  "tf"]


core_tfs
#%%

collection_dict = {}
for der_tf in der_tfs:
    for core_tf in core_tfs:

        pair = core_tf + "-" + der_tf

        results = test_cooccurrence(core_tf, derived_tf, cmplx)

        collection_dict[pair] = results
