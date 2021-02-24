import matplotlib.pyplot as plt
import numpy as np
import os, sys
import pandas as pd
import seaborn as sns
from scipy import stats


MPRAPATH = "/dors/capra_lab/projects/enhancer_ages/ernst16/new_data"

RE = "/dors/capra_lab/projects/enhancer_ages/ernst16/new_results/"

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
    sns.barplot(x = x, y = y, data = data,  order = order, n_boot=5000
    #showfliers = False, notch = True
    )

    counts = df.groupby("arch")[y].count()
    means = df.groupby("arch")[y].mean()
    derived_activity = df.loc[df.arch == "complex_derived", y]
    complex_activity = df.loc[df.arch == "complex_core", y]

    stat, p = stats.mannwhitneyu(derived_activity, complex_activity)
    xticklabs = ["simple", "complex\ncore", "complex\nderived"]

    ax.set(
    xticklabels = xticklabs,
    xlabel = "core v. der mwu p = %s\n%s" % (p, means),
    ylabel = "predicted activity (bp)",
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


def plot_dist(active_df):
    fig, ax = plt.subplots()
    for arch in active_df.arch.unique():
        test = active_df.loc[active_df.arch == arch, 'activity']
        sns.kdeplot(test, label = arch)
    ax.legend()

#%%

CELL_MODEL = "K562"
MPRAFILE = CELL_MODEL + "_combined_base_tiles_ages.bed"

MPRAF = os.path.join(MPRAPATH, MPRAFILE)

k562_info, k562_active_df, k562_wtact = formatdf(MPRAF, CELL_MODEL)

#%%

plot_dist(k562_active_df)

#%% plot results
CELL_MODEL = "K562"
x = "arch"
y = "act/synlen"
outf = make_pdf("%s_ernst_active_bases_dist.pdf"% CELL_MODEL, RE)
plot_activity(x, y, k562_active_df, outf, CELL_MODEL)


#%% weighted activity

CELL_MODEL = "K562"
x = "arch"
y = "max_act/synlen"
outf = make_pdf("%s_ernst_weighted_activity_bases_dist.pdf"% CELL_MODEL, RE)
plot_activity(x, y, k562_wtact, outf, CELL_MODEL)


#%% how many arch bases are there?

k562_freq = get_act_freq(k562_info, k562_active_df)
k562_freq
# about 2.5% of basepairs per architecture have activity >=2

#%% ### HEPG2 ###

CELL_MODEL = "HEPG2"
MPRAFILE = CELL_MODEL + "_combined_base_tiles_ages.bed"

MPRAF = os.path.join(MPRAPATH, MPRAFILE)

hepg2_info, hepg2_active_df, wtact_hepg2 = formatdf(MPRAF, CELL_MODEL)


#%%
plot_dist(hepg2_active_df)

#%% plot results
CELL_MODEL = "HEPG2"
x = "arch"
y = "act/synlen"
outf = make_pdf("%s_ernst_active_bases_dist.pdf"% CELL_MODEL, RE)
plot_activity(x, y, hepg2_active_df, outf, CELL_MODEL)
hepg2_active_df.groupby("arch")["activity"].mean()


#%%
CELL_MODEL = "HEPG2"
x = "arch"
y = "max_act/synlen"
outf = make_pdf("%s_ernst_weighted_activity_bases_dist.pdf"% CELL_MODEL, RE)
plot_activity(x, y, wtact_hepg2, outf, CELL_MODEL)

#%% how many arch bases are there?

hepg2_freq = get_act_freq(hepg2_info, hepg2_active_df)
hepg2_freq
# about 2.5% of basepairs per architecture have activity >=2
