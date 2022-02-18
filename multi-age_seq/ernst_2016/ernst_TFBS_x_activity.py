import os, sys
import pandas as pd
import matplotlib.pyplot as plt

MPRAPATH = "/dors/capra_lab/projects/enhancer_ages/ernst16/new_data"
CELL_MODEL = "HEPG2"
MPRAFILE = CELL_MODEL + "_combined_base_tiles_ages.bed"
MPRAF = os.path.join(MPRAPATH, MPRAFILE)

RE = "/dors/capra_lab/projects/enhancer_ages/landscape/results/ernst16/TFBS/"
if os.path.exists(RE) == False:
    os.mkdir(RE)
colors = [ "amber", "dusty purple", "windows blue"]
palette = sns.xkcd_palette(colors)
sns.palplot(palette)


#%% FUNCTIONS


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

    # normalize activity over syntenic length.
    cldf["act/synlen"] = cldf.activity.divide(cldf.syn_len)

    # limit analysis to bases with activity scores greater than or equal to 1
    active_cldf = cldf.loc[cldf.activity >= 1].copy()


    wtact = cldf.groupby(["syn_id", "syn_len", "arch"])["activity"].max().reset_index()
    wtact["max_act/synlen"] = wtact.activity.divide(wtact.syn_len)

    return cldf, active_cldf, wtact


def tfbs_chip_intersection(activity_f, cell_model, path):

    if cell_model == "HEPG2":
        cell_model = "HepG2"
    # save results to a file
    outfile = f"{path}/{cell_model}_ernst_x_ENCODE3_TFBS.bed"

    # get the matching encode cell line TFBS data.
    ENCODEPATH = "/dors/capra_lab/data/encode/encode3_hg38/TF/liftOver_hg19/cells/"
    encodefile = f"{cell_model}.bed.gz"
    encodef = os.path.join(ENCODEPATH, encodefile)

    # do bedtools intersection and write all overlapping/non-overlapping loci
    cmd = f"bedtools intersect -a {activity_f} -b {encodef} -wao > {outfile}"
    subprocess.call(cmd, shell = True)

    print(cmd)
    return outfile

def format_intersectiondf(int_f, info_df):

    # column name formatting...
    column_names = list(info_df) # get column names
    column_names.append("tf") # add TF column name
    usecols1 = [i for i in range(22)] # columns to select
    usecols2 = [27] # more columns to select (from TFBS intersection)
    usecols = usecols1 + usecols2 # combine the columns to select

    # let's open the file!
    df = pd.read_csv(int_f,
    header = None,
    usecols = usecols, sep ='\t', names = column_names)

    return df



def get_xlabels(data_all, data_act):

    ax1= data_all.groupby("arch")["enh_id"].count().reset_index().sort_values(by = "arch")
    ax2= data_act.groupby("arch")["enh_id"].count().reset_index().sort_values(by = "arch")

    archs = ["simple", "complex_core", "complex_derived"]

    ax1_lab = []
    ax2_lab = []

    for arch in archs:
        if arch in ax1.arch.to_list():
            n1 = ax1.loc[ax1.arch == arch, "enh_id"].iloc[0]
        else:
            n1 = 0
        x1_lab = f"{arch}\nn = {n1}"
        ax1_lab.append(x1_lab)


        if arch in ax2.arch.to_list():
            n2 = ax2.loc[ax2.arch == arch, "enh_id"].iloc[0]
        else:
            n2 = 0
        x2_lab = f"{arch}\nn = {n2}"
        ax2_lab.append(x2_lab)

    return ax1_lab, ax2_lab


def plot_tf_activity(tf, intdf, RE):


    tfdf = intdf.loc[intdf.tf == tf] # get the TF-specific dataframe

    x, y = "arch", "act/synlen"

    data_all = tfdf
    data_act = tfdf.loc[tfdf.activity >=1]
    ax1_lab, ax2_lab = get_xlabels(data_all, data_act)
    ax1_lab
    order =["simple", "complex_core", "complex_derived"]
    fig, (ax1, ax2) = plt.subplots(ncols = 2, figsize = (16, 6))

    sns.barplot(
    x=x,
    y=y,
    data = data_all,
    ax = ax1,
    order = order,
    palette = palette)

    ax1.set(title = f"all_enh-{tf}", xlabel = '')
    ax1.set_xticklabels(ax1_lab, rotation = 90)

    sns.barplot(
    x=x,
    y=y,
    data = data_act,
    ax = ax2,
    order = order,
    palette = palette)

    ax2.set(title = f"active_enh-{tf}",  xlabel = '')
    ax2.set_xticklabels(ax2_lab, rotation = 90)

    outf = os.path.join(RE, f"{tf}.pdf")
    plt.savefig(outf, bbox_inches = "tight")

#%% RUN ANALYSIS


info_df, active_df, wt_act = formatdf(MPRAF, CELL_MODEL)

# prepare ernst file to intersect w/ ENCODE TFBS file

activity_f = os.path.join(MPRAPATH, f"{CELL_MODEL}_activity_df.bed")

info_df.to_csv(activity_f, sep = '\t', header = None, index = False) # save the dataframe for intersection w/ TFBS

int_f = tfbs_chip_intersection(activity_f, CELL_MODEL, MPRAPATH) # intersect with TFBS

#%% GET intersection df


intdf = format_intersectiondf(int_f, info_df)

intdf.head()

#%%
tf_list = list(intdf.tf.unique())
for tf in tf_list[47:]:
    if tf !=".":
        print(tf)
        plot_tf_activity(tf, intdf, RE)
