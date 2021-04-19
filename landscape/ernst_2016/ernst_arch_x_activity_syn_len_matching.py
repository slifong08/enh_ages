import matplotlib.pyplot as plt
import numpy as np
import os, sys
import pandas as pd
import seaborn as sns
from scipy import stats


MPRAPATH = "/dors/capra_lab/projects/enhancer_ages/ernst16/new_data"

RE = "/dors/capra_lab/projects/enhancer_ages/landscape/results/ernst16/"

if os.path.exists(RE)== False:
    os.mkdir(RE)

colors = [ "amber", "dusty purple", "windows blue"]
palette = sns.xkcd_palette(colors)
sns.palplot(palette)
#%%


def formatdf(f, cell_model):

    if cell_model == "HEPG2":
        cell_model = "Hepg2"

    # assign column names
    cols = ["chr_bp", "start_bp", "end_bp", "activity",
    "chHMM_state", "tile_id",
    "chr_syn", "start_syn", "end_syn",
    "enh_id", "chr_enh", "start_enh", "end_enh",
    "seg_index", "core_remodeling", "core", "mrca"]

    # open dataframe
    df = pd.read_csv(f, sep = '\t', header = None, names = cols)

    df["cell_enh"] = df.tile_id.apply(lambda x: x.split("_")[0])


    # keep only tiles designed from cis-cell_model, only enhancer ChromHMM states
    cldf = df.loc[(df.cell_enh.str.contains(cell_model)) &
    (df.chHMM_state>=5) & (df.chHMM_state<=8)].copy()

    # assign architecture
    cldf["arch"] = "complex_core"
    cldf.loc[cldf.core == 0, "arch"] = "complex_derived"
    cldf.loc[cldf.core_remodeling == 0, "arch"] = "core"

    cldf["syn_id"] = cldf.chr_syn + ":" +cldf.start_syn.map(str) + "-" + cldf.end_syn.map(str)
    cldf["syn_len"] = cldf.end_syn - cldf.start_syn

    # add age categories
    syn_gen_bkgd_file = "/dors/capra_lab/projects/enhancer_ages/hg19_syn_gen_bkgd.tsv"
    syn_gen_bkgd= pd.read_csv(syn_gen_bkgd_file, sep = '\t') # read the file
    syn_gen_bkgd[["mrca", "mrca_2"]] = syn_gen_bkgd[["mrca", "mrca_2"]].round(3) # round the ages

    syn_gen_bkgd = syn_gen_bkgd[["mrca", "taxon", "mrca_2", "taxon2"]] # whittle down the df

    # round the ages
    cldf["mrca"] = cldf["mrca"].round(3)
    cldf = pd.merge(cldf, syn_gen_bkgd, how = "left", on = "mrca")
    cldf = cldf.drop(["mrca"], axis = 1).drop_duplicates()

    core_age = cldf.groupby("enh_id")["mrca_2"].max().reset_index()
    core_age.columns = ["enh_id", "core_mrca_2"]
    cldf = pd.merge(cldf, core_age, how = "left", on = "enh_id")

    # limit analysis to bases with activity
    active_cldf = cldf.loc[cldf.activity >= 1].copy()
    active_cldf = cldf.loc[(cldf.activity >= 1) & (cldf.syn_len > 50)].copy() # require that syntenic blocks are at least 50bp long.
    active_cldf["act/synlen"] = active_cldf.activity.divide(active_cldf.syn_len)

    wtact = cldf.groupby(["syn_id", "syn_len", "arch"])["activity"].max().reset_index()
    wtact["max_act/synlen"] = wtact.activity.divide(wtact.syn_len)

    return cldf, active_cldf, wtact


def custom_round(x, base=10):
    return int(base * round(float(x)/base))


def match_len(core_df, der_df, base_len):

    columns = ["syn_id", "syn_len"]
    columns_names = ["matching_ids", "matching_len"]

    core = core_df[columns].drop_duplicates()

    core.columns = columns_names
    core.matching_len = core.matching_len.astype(float).apply(lambda x: custom_round(x, base=base_len)) # round to the nearest 100bp

    der = der_df[columns].drop_duplicates()
    der.columns = columns_names

    der.matching_len = der.matching_len.astype(float).apply(lambda x: custom_round(x, base=base_len))

    lens = set(list(core.matching_len.unique()) + list(der.matching_len.unique())) # find the intersecting lengths

    match_dict = {}

    for length in lens:
        dern = der.loc[der.matching_len == length].size
        coren = core.loc[core.matching_len == length].size

        sn = min(dern, coren)

        if length > 0 and sn > 0:


            # find 100 matched enhancer samples

            der_ids = der.loc[der.matching_len == length].sample(n = sn, replace = True) # sample w/ replacement
            core_ids = core.loc[core.matching_len == length].sample(n = sn, replace = True) # sample w/ replacement
            balanced = pd.concat([core_ids, der_ids])
            match_dict[sn] = balanced

    final_matched_id = pd.concat(match_dict.values())

    return final_matched_id.matching_ids.unique()


def plot_activity(x, y, df, outf, cell_model):


    order = ["core", "complex_core", "complex_derived"]
    data = df

    fig, ax = plt.subplots(figsize = (6,6))
    sns.set("poster")
    sns.set_style("white")
    sns.barplot(x = x, y = y, data = data,  order = order, palette = palette)
    #showfliers = False, notch = True


    counts = df.groupby("arch")[y].count()
    means = df.groupby("arch")[y].mean()
    derived_activity = df.loc[df.arch == "complex_derived", y]
    complex_activity = df.loc[df.arch == "complex_core", y]

    stat, p = stats.mannwhitneyu(derived_activity, complex_activity)
    xticklabs = ["core", "der\ncore", "der\nderived"]

    ax.set(
    xticklabels = xticklabs,
    xlabel = "core v. der mwu p = %s\n%s" % (p, means),
    ylabel = "predicted activity (bp)/syn len",
    title = cell_model
    )


    plt.savefig(outf, bbox_inches = "tight")


def get_act_freq(info_df, active_df):

    n_total = info_df.groupby("arch")["enh_id"].count().reset_index()
    n_total.columns = ["arch", "total_bp_count"]
    n_active = active_df.groupby("arch")["enh_id"].count().reset_index()
    n_active.columns = ["arch", "active_bp_count"]

    freq = pd.merge(n_total, n_active)
    freq["freq"] = freq.active_bp_count.divide(freq.total_bp_count)

    return freq


def plot_dist(active_df, variable):
    fig, ax = plt.subplots()
    for arch in active_df.arch.unique():
        test = active_df.loc[active_df.arch == arch, variable]
        sns.kdeplot(test, label = arch)
    ax.legend()


def plot_stratified_age(act_df, outf, x):


    y = "act/synlen"
    hue = "arch"
    data = act_df
    order = ["core", "complex_core", "complex_derived"]


    fig, ax = plt.subplots(figsize = (6,6))
    sns.set("poster")
    sns.barplot(x = x,
    y = y,
    data = data,
    hue = hue,
    hue_order = order,
    palette = palette
    )
    labs = ["prim", "euar", "bore", "euth", "ther", "mam", "amni", "tetr", "vert"]
    ax.set_xticklabels(labs, rotation = 90)
    ax.legend(bbox_to_anchor = (1,1))

    plt.savefig(outf, bbox_inches = "tight")


#%%

CELL_MODEL = "K562"
MPRAFILE = CELL_MODEL + "_combined_base_tiles_ages.bed"

MPRAF = os.path.join(MPRAPATH, MPRAFILE)

k562_info, k562_active_df, k562_wtact = formatdf(MPRAF, CELL_MODEL)

#%%

core_df = k562_active_df.loc[k562_active_df.arch == "complex_core"]
der_df = k562_active_df.loc[k562_active_df.arch == "complex_derived"]

base_len = 10
matched_ids = match_len(core_df, der_df, base_len)

matched_k562 = k562_active_df.loc[k562_active_df["syn_id"].isin(matched_ids)]
#%%
sns.distplot(k562_active_df.loc[k562_active_df.arch == "complex_core", "syn_len"], label = "core")
sns.distplot(k562_active_df.loc[k562_active_df.arch == "complex_derived", "syn_len"], label = "der")

plt.legend()
outf = f"{RE}{CELL_MODEL}_ernst_active_bases_dist_NOT_len_matched.pdf"
plt.savefig(outf, bbox_inches = "tight")

#%% plot activity stratified by age
sns.distplot(matched_k562.loc[matched_k562.arch == "complex_core", "syn_len"], label = "core")
sns.distplot(matched_k562.loc[matched_k562.arch == "complex_derived", "syn_len"], label = "der")
plt.legend()
outf =  f"{RE}{CELL_MODEL}_ernst_active_bases_dist_len_matched.pdf"
plt.savefig(outf, bbox_inches = "tight")
#%%
CELL_MODEL = "K562"

outf =f"{RE}{CELL_MODEL}_ernst_active_bases_dist_mrca2_len_matched.pdf"
x = "mrca_2"
plot_stratified_age(matched_k562, outf, x)


outf = f"{RE}{CELL_MODEL}_ernst_active_bases_dist_core_mrca2_len_matched.pdf"
x = "core_mrca_2"
plot_stratified_age(matched_k562, outf, x)


#%%

plot_dist(matched_k562, "activity")

plot_dist(matched_k562, "act/synlen")

#%% plot results

CELL_MODEL = "K562"
x = "arch"
y = "act/synlen"
outf = f"{RE}{CELL_MODEL}_ernst_active_bases_dist_mrca2_len_matched.pdf"
plot_activity(x, y, matched_k562, outf, CELL_MODEL)
matched_k562.groupby("arch")[y].mean()

"""
arch            act/synlen
complex_core       0.009983
complex_derived    0.014807
core             0.006812
"""



#%% how many arch bases are there?

k562_freq = get_act_freq(k562_info, k562_active_df)
k562_freq
# about 2.5% of basepairs per architecture have activity >=2 (data not shown)
# about 7.5 - 8% of basepairs per architecture have activity >=1 (data not shown)
"""
	arch	total_bp_count	active_bp_count	freq
0	complex_core	123074	9906	0.080488
1	complex_derived	108501	8119	0.074829
2	core	343675	26658	0.077567
"""
#%% ### HEPG2 ###

CELL_MODEL = "HEPG2"
MPRAFILE = CELL_MODEL + "_combined_base_tiles_ages.bed"

MPRAF = os.path.join(MPRAPATH, MPRAFILE)

hepg2_info, hepg2_active_df, wtact_hepg2 = formatdf(MPRAF, CELL_MODEL)


#%%

core_df = hepg2_active_df.loc[hepg2_active_df.arch == "complex_core"]
der_df = hepg2_active_df.loc[hepg2_active_df.arch == "complex_derived"]

base_len = 10
matched_ids = match_len(core_df, der_df, base_len)

matched_hepg2 = hepg2_active_df.loc[hepg2_active_df["syn_id"].isin(matched_ids)]
#%%
x = "mrca_2"
outf = f"{RE}{CELL_MODEL}_ernst_active_bases_dist_mrca2_len_matched.pdf"
plot_stratified_age(matched_hepg2, outf, x)

x = "core_mrca_2"
outf = f"{RE}{CELL_MODEL}_ernst_active_bases_dist_core_mrca2_len_matched.pdf"
plot_stratified_age(matched_hepg2, outf, x)

#%%
plot_dist(matched_hepg2, "activity")
plot_dist(matched_hepg2, "act/synlen")

#%% plot results
CELL_MODEL = "HEPG2"
x = "arch"
y = "act/synlen"
outf = f"{RE}{CELL_MODEL}_ernst_active_bases_dist_len_matched.pdf"
plot_activity(x, y, matched_hepg2, outf, CELL_MODEL)
matched_hepg2.groupby("arch")[y].mean()

"""
arch        activity/syn_len
complex_core       0.010030
complex_derived    0.013797
core             0.006812
"""
#%%
sns.distplot(hepg2_active_df.loc[hepg2_active_df.arch == "complex_core", "syn_len"], label = "core")
sns.distplot(hepg2_active_df.loc[hepg2_active_df.arch == "complex_derived", "syn_len"], label = "der")

plt.legend()
outf = f"{RE}{CELL_MODEL}_ernst_active_bases_dist_NOT_len_matched.pdf"
plt.savefig(outf, bbox_inches = "tight")


#%% plot activity stratified by age
sns.distplot(matched_hepg2.loc[matched_hepg2.arch == "complex_core", "syn_len"], label = "core")
sns.distplot(matched_hepg2.loc[matched_hepg2.arch == "complex_derived", "syn_len"], label = "der")
plt.legend()
outf = f"{RE}{CELL_MODEL}_ernst_active_bases_dist_len_matched.pdf"
plt.savefig(outf, bbox_inches = "tight")
#%% how many arch bases are there?

hepg2_freq = get_act_freq(hepg2_info, hepg2_active_df)
hepg2_freq
# about 2.5% of basepairs per architecture have activity >=2 (data not shown)
# about 6.6 - 8% of basepairs per architecture have activity >=1
"""
	arch	total_bp_count	active_bp_count	freq
0	complex_core	133777	9076	0.067844
1	complex_derived	103403	6498	0.062842
2	core	345740	27672	0.080037
"""
#%%

hepg2_active_df.head()
cols = ["enh_id", "core", "arch", "syn_len", "mrca_2", "taxon2", "core_mrca_2", "act/synlen"]
cores = hepg2_active_df.loc[hepg2_active_df.arch == "complex_core", cols].drop_duplicates()
ders = hepg2_active_df.loc[hepg2_active_df.arch == "complex_derived", cols].drop_duplicates()

arch = pd.merge(cores, ders, how = "outer", on = "enh_id")
arch.shape
arch = arch.dropna()
arch.shape
arch.head()
x = "act/synlen_x"
y = "act/synlen_y"
data = arch
sns.jointplot(x = x, y = y, data = data)
#%%

from sklearn.linear_model import LinearRegression

# does core activity predict derived activity?
X = np.array(arch["act/synlen_x"]).reshape(-1, 1)
y = np.array(arch["act/synlen_y"]).reshape(-1, 1)
reg = LinearRegression().fit(X, y)

reg.score(X, y) # 0.00963961057383722 kind of

reg.coef_ #0.13432273

# does derived activity predict core activity?
# does core activity predict derived activity?
y = np.array(arch["act/synlen_x"]).reshape(-1, 1)
x = np.array(arch["act/synlen_y"]).reshape(-1, 1)
x
reg = LinearRegression().fit(X, y)

reg.score(X, y) # 0.00963961057383722 kind of

reg.coef_ #0.13432273
