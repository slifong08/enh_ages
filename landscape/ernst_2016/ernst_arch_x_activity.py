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
    cldf.loc[cldf.core_remodeling == 0, "arch"] = "simple"

    cldf["syn_id"] = cldf.chr_syn + ":" +cldf.start_syn.map(str) + "-" + cldf.end_syn.map(str)
    cldf["syn_len"] = cldf.end_syn - cldf.start_syn

    # limit analysis to bases with activity
    active_cldf = cldf.loc[cldf.activity >= 1].copy()
    active_cldf["act/synlen"] = active_cldf.activity.divide(active_cldf.syn_len)

    wtact = cldf.groupby(["syn_id", "syn_len", "arch"])["activity"].max().reset_index()
    wtact["max_act/synlen"] = wtact.activity.divide(wtact.syn_len)

    return cldf, active_cldf, wtact


def make_pdf(file_name, RE):

    OUTFILE = file_name + ".pdf"
    OUTF = os.path.join(RE, OUTFILE)

    return OUTF


def plot_activity(x, y, df, outf, cell_model):


    order = ["simple", "complex_core", "complex_derived"]
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
    xticklabs = ["simple", "complex\ncore", "complex\nderived"]

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

#%%

CELL_MODEL = "K562"
MPRAFILE = CELL_MODEL + "_combined_base_tiles_ages.bed"

MPRAF = os.path.join(MPRAPATH, MPRAFILE)

k562_info, k562_active_df, k562_wtact = formatdf(MPRAF, CELL_MODEL)

#%%

plot_dist(k562_active_df, "activity")

plot_dist(k562_active_df, "act/synlen")

#%% plot results

CELL_MODEL = "K562"
x = "arch"
y = "act/synlen"
outf = make_pdf("%s_ernst_active_bases_dist.pdf"% CELL_MODEL, RE)
plot_activity(x, y, k562_active_df, outf, CELL_MODEL)
k562_active_df.groupby("arch")[y].mean()

"""
arch            act/synlen
complex_core       0.009983
complex_derived    0.014807
simple             0.006812
"""

#%% max weighted activity

CELL_MODEL = "K562"
x = "arch"
y = "max_act/synlen"
outf = make_pdf("%s_ernst_max_activity_per_bases_dist.pdf"% CELL_MODEL, RE)
plot_activity(x, y, k562_wtact, outf, CELL_MODEL)


#%% how many arch bases are there?

k562_freq = get_act_freq(k562_info, k562_active_df)
k562_freq
# about 2.5% of basepairs per architecture have activity >=2 (data not shown)
# about 7.5 - 8% of basepairs per architecture have activity >=1 (data not shown)
"""
	arch	total_bp_count	active_bp_count	freq
0	complex_core	123074	9906	0.080488
1	complex_derived	108501	8119	0.074829
2	simple	343675	26658	0.077567
"""
#%% ### HEPG2 ###

CELL_MODEL = "HEPG2"
MPRAFILE = CELL_MODEL + "_combined_base_tiles_ages.bed"

MPRAF = os.path.join(MPRAPATH, MPRAFILE)

hepg2_info, hepg2_active_df, wtact_hepg2 = formatdf(MPRAF, CELL_MODEL)


#%%
plot_dist(hepg2_active_df, "activity")
plot_dist(hepg2_active_df, "act/synlen")

#%% plot results
CELL_MODEL = "HEPG2"
x = "arch"
y = "act/synlen"
outf = make_pdf("%s_ernst_active_bases_dist.pdf"% CELL_MODEL, RE)
plot_activity(x, y, hepg2_active_df, outf, CELL_MODEL)
hepg2_active_df.groupby("arch")[y].mean()

"""
arch        activity/syn_len
complex_core       0.010030
complex_derived    0.013797
simple             0.006812
"""
#%%
CELL_MODEL = "HEPG2"
x = "arch"
y = "max_act/synlen"
outf = make_pdf("%s_ernst_weighted_activity_bases_dist.pdf"% CELL_MODEL, RE)
plot_activity(x, y, wtact_hepg2, outf, CELL_MODEL)

#%% how many arch bases are there?

hepg2_freq = get_act_freq(hepg2_info, hepg2_active_df)
hepg2_freq
# about 2.5% of basepairs per architecture have activity >=2 (data not shown)
# about 6.6 - 8% of basepairs per architecture have activity >=1
"""
	arch	total_bp_count	active_bp_count	freq
0	complex_core	133777	9090	0.067949
1	complex_derived	103403	6894	0.066671
2	simple	345740	27672	0.080037
"""
